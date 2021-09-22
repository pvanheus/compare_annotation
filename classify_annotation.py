#!/usr/bin/env python

from __future__ import annotations
import argparse
import dataclasses
from itertools import combinations
import sys
from typing import TextIO

from intervaltree import IntervalTree, Interval
import CPT_GFFParser


@dataclasses.dataclass
class IntervalData:
    name: str
    sub_intervals: IntervalTree


def extract_interval(tree: IntervalTree, interval: Interval) -> Interval:
    """ Find the exact interval matching coordinates of 'interval' """
    interval_list = [i for i in tree.envelop(interval) if i.begin == interval.begin and i.end == interval.end]
    return interval_list[0]


def contains_interval(tree: IntervalTree, interval: Interval) -> bool:
    """Test if an interval tree contains an interval

    Arguments:
    tree -- the IntervalTree being searched
    interval -- the Interval being searched for in the tree. Only begin and end is considered, not the data part of the Interval
    """
    intervals = tree.envelop(interval)
    for known_interval in intervals:
        if known_interval.begin == interval.begin and known_interval.end == interval.end:
            return True
    else:
        return False


def gff_to_intervals(gff_file: TextIO, use_exon: bool = True) -> tuple[IntervalTree, IntervalTree]:
    """ Convert a GFF3 file into two IntervalTrees with genes and their exons """
    gene_tree_forward = IntervalTree()
    gene_tree_reverse = IntervalTree()

    parser = CPT_GFFParser.gffParse(gff_file)
    start = end = -1
    for record in parser:
        for feature in record.features:
            if feature.type != "gene" or len(feature.sub_features) > 0 and feature.sub_features[0].type != "mRNA":
                continue
            gene_start = feature.location.start
            gene_end = feature.location.end
            strand = feature.location.strand
            exons = IntervalData(feature.qualifiers.get('Name', feature.id)[0], IntervalTree())
            gene_interval = Interval(gene_start, gene_end, exons)
            for subfeature in feature.sub_features:
                if subfeature.type != 'mRNA':
                    continue
                # we're now looking at transcripts... we want to capture each exon of this gene
                # mrna_interval = Interval(subfeature.location.start, subfeature.location.end, exons)
                transcript_name = subfeature.id
                ssf_list = [(ssf.location.start, ssf.location.end, ssf.type, ssf.id, ssf.strand)
                            for ssf in subfeature.sub_features]
                ssf_list_len = len(ssf_list)
                for i in range(ssf_list_len):
                    (start, end, feature_type, id, strand) = ssf_list[i]
                    save_feature = False
                    if use_exon and feature_type == 'exon':
                        save_feature = True
                    elif not use_exon and feature_type == 'CDS':
                        save_feature = True
                        if i > 0 and ssf_list[i-1][2] == 'five_prime_UTR':
                            if strand == 1:
                                start = ssf_list[i-1][0]
                            elif strand == -1:
                                end = ssf_list[i-1][1]
                        if i < (ssf_list_len - 1) and ssf_list[i+1][2] == 'three_prime_UTR':
                            if strand == 1:
                                end = ssf_list[i+1][1]
                            elif strand == -1:
                                start = ssf_list[i+1][0]
                    if save_feature:
                        # save the exon to the exon list, keep track of its associated transcript ID
                        exon_interval = Interval(start, end, [transcript_name])
                        if contains_interval(exons.sub_intervals, exon_interval):
                            existing_interval = extract_interval(exons.sub_intervals, exon_interval)
                            existing_interval.data.append(transcript_name)
                        else:
                            exons.sub_intervals.add(exon_interval)
            if strand == 1:
                gene_tree_forward.add(gene_interval)
            elif strand == -1:
                gene_tree_reverse.add(gene_interval)
            else:
                raise ValueError("Unknown strand: " + str(strand))

    return (gene_tree_forward, gene_tree_reverse)


def flatten_set(nested_set: set[tuple[str, str]]) -> set[str]:
    """ Given a set of tuples, return a set with the members of the tuples """
    flat_set = set()
    for pair in nested_set:
        for element in pair:
            flat_set.add(element)
    return flat_set


def classify_annotation(tree):
    for gene_interval in tree:
        # compute M-N-O
        m_num = n_num = o_num = 0
        # transcript pairs that share exons - the number of distinct pairs is O in the M-N-O scheme
        transcripts_with_same_exons = set()
        transcripts_with_overlapping_exons = set()
        all_transcript_names = set()
        for exon_interval in gene_interval.data.sub_intervals:
            for transcript_name in exon_interval.data:
                all_transcript_names.add(transcript_name)
            if len(exon_interval.data) > 1:
                transcript_names = sorted(exon_interval.data)
                for combination in combinations(transcript_names, 2):
                    transcripts_with_same_exons.add(combination)

            overlaps = gene_interval.data.sub_intervals.overlap(exon_interval)
            if len(overlaps) > 1:
                # because of how intervals are captures this only counts overlapping but not identical exons
                # the number of transcripts with overlapping exons is N in the M-N-O scheme
                overlaps.remove(exon_interval)
                # overlaps now contains only non-identical exons
                transcript_names = set([exon_interval.data[0]])
                for interval in overlaps:
                    # add all the transcripts that have this exon - there could be more than one
                    for transcript_name in interval.data:
                        transcript_names.add(transcript_name)
                for combination in combinations(sorted(transcript_names), 2):
                    transcripts_with_overlapping_exons.add(combination)

        transcript_pairs = set()  # these are pairs of transcripts with no overlapping or common exons
        if len(all_transcript_names) > 1:
            for combination in combinations(sorted(all_transcript_names), 2):
                transcript_pairs.add(combination)
            for pair in transcripts_with_same_exons:
                if pair in transcripts_with_overlapping_exons:
                    transcripts_with_overlapping_exons.remove(pair)
                if pair in transcript_pairs:
                    transcript_pairs.remove(pair)
            for pair in transcripts_with_overlapping_exons:
                if pair in transcript_pairs:
                    transcript_pairs.remove(pair)

            m_num = len(transcript_pairs)
            n_num = len(transcripts_with_overlapping_exons)
            o_num = len(transcripts_with_same_exons)
        print(gene_interval.data.name, len(gene_interval.data.sub_intervals), m_num, n_num, o_num)


parser = argparse.ArgumentParser()
parser.add_argument('--use_exon', default=False, action='store_true')
parser.add_argument('gff_file', type=argparse.FileType())
args = parser.parse_args()

if not args.use_exon:
    print("Warning: using CDS rather than exon annotation has not been well tested", file=sys.stderr)

(fwd_tree, rev_tree) = gff_to_intervals(args.gff_file, args.use_exon)
print("Forward")
classify_annotation(fwd_tree)
print("Reverse")
classify_annotation(rev_tree)
