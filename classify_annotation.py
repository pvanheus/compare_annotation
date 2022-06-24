#!/usr/bin/env python

from __future__ import annotations
import argparse
import dataclasses
from itertools import combinations
import sys
import textwrap
from typing import TextIO, Union

from intervaltree import IntervalTree, Interval
import CPT_GFFParser


@dataclasses.dataclass
class IntervalData:
    name: str
    sub_intervals: Union[IntervalTree, list[IntervalTree]]


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


get_start = lambda x: x[0]
get_end = lambda x: x[1]
get_type = lambda x: x[2]

# TODO: write tests!
def gff_to_intervals(gff_file: TextIO, use_exon: bool = True, add_transcripts = False) -> tuple[IntervalTree, IntervalTree]:
    """ Convert a GFF3 file into two IntervalTrees with genes and their exons """
    five_prime_utr_str = 'five_prime_UTR'
    three_prime_utr_str = 'three_prime_UTR'

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
            exons = IntervalData(feature.qualifiers.get('ID', feature.id)[0], IntervalTree())
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
                # five_prime_utr_count = 0
                for i in range(ssf_list_len):
                    (start, end, feature_type, id, strand) = ssf_list[i]
                    save_feature = False
                    if use_exon and feature_type == 'exon':
                        save_feature = True
                    elif not use_exon and (feature_type == 'CDS' or feature_type == five_prime_utr_str or feature_type == three_prime_utr_str):
                        save_feature = True
                        if feature_type == five_prime_utr_str and get_type(ssf_list[i+1]) == 'CDS':
                            # don't save a five_prime_UTR followed by a CDS - the 5' UTR will be used to extend the CDS
                            save_feature = False
                        elif feature_type == three_prime_utr_str and get_type(ssf_list[i-1]) == 'CDS':
                            # don't save a three_prime_UTR that follows a CDS - the 3' UTR will be used to extend the CDS
                            save_feature = False
                        # If we are using a form of GFF3 that does not annotate exons,
                        # we have to use CDS, three_prime and five_prime UTR features.
                        # print(strand, transcript_name, feature_type, save_feature, i, ssf_list_len - 1)
                        # if ssf_list[i][2] == 'five_prime_UTR':
                        #     five_prime_utr_count += 1
                        feature_type = get_type(ssf_list[i])
                        if i > 0 and feature_type == 'CDS' and get_type(ssf_list[i-1]) == five_prime_utr_str:
                            # five_prime_UTR followed by exon - extend exon start "backwards"
                            if strand == 1:
                                start = get_start(ssf_list[i-1])
                            elif strand == -1:
                                end = get_end(ssf_list[i-1])
                        if i < (ssf_list_len - 1) and feature_type == 'CDS' and get_type(ssf_list[i+1]) == three_prime_utr_str:
                            # exon followed by three_prime_UTR - extend exon end forwards
                            if strand == 1:
                                end = get_end(ssf_list[i+1])
                            elif strand == -1:
                                start = get_start(ssf_list[i+1])
                        # elif strand == -1:
                        #     if i < (ssf_list_len - 1) and feature_type == 'CDS' and ssf_list[i+1][2] == 'five_prime_UTR':
                        #         # exon followed by five_prime_UTR - extend exon end forwards
                        #         end = ssf_list[i+1][1]
                        #     elif i > 0 and feature_type == 'CDS' and ssf_list[i-1][2] == 'three_prime_UTR':
                        #         start = ssf_list[i-1][0]
                    if save_feature:
                        # save the exon to the exon list, keep track of its associated transcript ID
                        exon_interval = Interval(start, end, [transcript_name])
                        if contains_interval(exons.sub_intervals, exon_interval):
                            existing_interval = extract_interval(exons.sub_intervals, exon_interval)
                            existing_interval.data.append(transcript_name)
                        else:
                            exons.sub_intervals.add(exon_interval)
                # if five_prime_utr_count > 0:
                #     print(f"Warning: {five_prime_utr_count} five_prime_UTR features found in {transcript_name}", file=sys.stderr)
            if strand == 1:
                gene_tree_forward.add(gene_interval)
            elif strand == -1:
                gene_tree_reverse.add(gene_interval)
            else:
                raise ValueError("Unknown strand: " + str(strand))

    if add_transcripts:
        gene_tree_forward = add_transcript_level(gene_tree_forward)
        gene_tree_reverse = add_transcript_level(gene_tree_reverse)
    return (gene_tree_forward, gene_tree_reverse)


def flatten_set(nested_set: set[tuple[str, str]]) -> set[str]:
    """ Given a set of tuples, return a set with the members of the tuples """
    flat_set = set()
    for pair in nested_set:
        for element in pair:
            flat_set.add(element)
    return flat_set


def classify_annotation(tree):
    """ classify annotation according to Eilbeck et al M-N-O scheme - doi:10.1186/1471-2105-10-67
        input is gene tree in gene -> exons form """
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
        if m_num == 0 and n_num == 0 and o_num == 0:
            # this is a singleton gene
            print(gene_interval.data.name, len(gene_interval.data.sub_intervals))
        else:
            print(gene_interval.data.name, len(gene_interval.data.sub_intervals), ':', m_num, n_num, o_num)


def add_transcript_level(gene_tree):
    """ take a gene tree with structure gene -> exons and add a transcript level, so its 
        gene -> list[transcript] -> exons"""
    for gene in gene_tree:
        # gene is an interval with exons in the gene_data
        exon_data = gene.data.sub_intervals
        transcripts = {}
        for exon in exon_data:
            for transcript_name in exon.data:
                if transcript_name not in transcripts:
                    transcripts[transcript_name] = IntervalTree()
                transcripts[transcript_name].add(exon)
        gene.data.sub_intervals = []
        for transcript_name in transcripts:
            transcript = IntervalData(transcript_name, transcripts[transcript_name])
            gene.data.sub_intervals.append(transcript)
    return gene_tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Classify gene annotation by splice complexity', 
                                     epilog=textwrap.dedent('''\
                                        For each gene, this script outputs either the number of exons, or, for genes
                                        with multiple isoforms, its classification according to the M-N-O scheme. 
                                        The M-N-O scheme is a way to classify the splice complexity of a gene. 
                                        M: the number of transcript pairs that don't share any exons,
                                        N: the number of transcript pairs that share overlapping exons but not exon boundaries and
                                        O: the number of transcript pairs that share at least some exons.'''))
    parser.add_argument('--use_exon', default=False, action='store_true', help='Use exon features instead of CDS')
    parser.add_argument('gff_file', type=argparse.FileType('r', encoding='latin-1'))
    args = parser.parse_args()

    if not args.use_exon:
        print("Warning: using CDS rather than exon annotation has not been well tested", file=sys.stderr)

    (fwd_tree, rev_tree) = gff_to_intervals(args.gff_file, args.use_exon)
    print("Forward")
    classify_annotation(fwd_tree)
    print()
    print("Reverse")
    classify_annotation(rev_tree)
