#!/usr/bin/env python

from __future__ import print_function, division
import argparse

from intervaltree import Interval, IntervalTree

from classify_annotation import gff_to_intervals, add_transcript_level

def compute_overlap(anno1: IntervalTree, anno2: IntervalTree) -> float:
    overall_length = 0
    overall_overlap = 0
    for interval in anno1:
        # note: intervals are "half open", i.e. begin is in the interval, end is outside
        # by contrast GFF3 intervals are open
        interval_length = interval.length()
        overall_length += interval_length
        # compute how much ove anno2 overlaps this interval in anno1
        overlap_distance = 0
        for overlap in anno2.overlap(interval):
            if interval.begin <= overlap.begin and overlap.end <= interval.end:
                # overlap is enveloped by interval
                overlap_distance += overlap.length()
            elif interval.begin <= overlap.begin and overlap.end > interval.end:
                # interval is before overlap
                overlap_distance += (interval.end - overlap.begin)
            elif overlap.begin < interval.begin and interval.begin < overlap.end:
                # interval is after overlap
                overlap_distance += (overlap.end - interval.begin)
            elif overlap.begin <= interval.begin and overlap.end <= interval.end:
                # interval is enclosed in overlap
                overlap_distance += interval_length
        overall_overlap += overlap_distance
    # overall_length is currently the length of all exons in transcript anno1
    # overall_overlap is the number of bases that anno1 overlaps anno2
    return overall_overlap / overall_length


def compute_discordance(anno1: IntervalTree, anno2: IntervalTree) -> float:
    # anno1 and anno2 are intervaltrees, each representing a transcript
    # the goal of this function is to report the discordance (D) between anno1 and anno2
    anno1_anno2_overlap = compute_overlap(anno1, anno2)
    anno2_anno1_overlap = compute_overlap(anno2, anno1)
    concordance = (anno1_anno2_overlap + anno2_anno1_overlap) / 2
    discordance = 1 - concordance
    return discordance


def compute_aed(anno1, anno2):
    """Compute the annotation edit distance between two sets of annotation.
       Each annotate set is a gene -> transcript -> exon tree as per those produced by gff_to_intervals.
       Genes are matched up in pairs and then AED is computed for each gene pair and returned.
       TODO: Find gene pairs using best BLAST hit, not overlapping annotation, as described in 
             Eilbeck 2009 'Tracking annotation from release to release'"""
    aed_by_gene = {}
    for gene_interval in anno1:
        for gene_interval2 in anno2:
            if gene_interval2.overlaps(gene_interval):
                transcript_pairs = []
                for transcript1 in gene_interval.data.sub_intervals:
                    min_discordance = 1
                    for transcript2 in gene_interval2.data.sub_intervals:
                        # each transcript is a IntervalTree representing a transcript
                        discordance = compute_discordance(transcript1.sub_intervals, transcript2.sub_intervals)
                        if discordance <= min_discordance:
                            pair = (transcript1, transcript2)
                            min_discordance = discordance
                    transcript_pairs.append(pair)
                AED = 0
                for pair in transcript_pairs:
                    AED += compute_discordance(pair[0].sub_intervals, pair[1].sub_intervals)
                AED = AED / len(transcript_pairs)
                aed_by_gene[gene_interval.data.name + '-' + gene_interval2.data.name] = AED
    return aed_by_gene


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--use_exon', default=False, action='store_true')
    parser.add_argument('anno1_file', type=argparse.FileType())
    parser.add_argument('anno2_file', type=argparse.FileType())
    args = parser.parse_args()

    (anno1_forward, anno1_reverse) = gff_to_intervals(args.anno1_file, args.use_exon, add_transcripts=True)
    (anno2_forward, anno2_reverse) = gff_to_intervals(args.anno2_file, args.use_exon, add_transcripts=True)
    
                # print(gene_interval2.data.sub_intervals)
                # print(gene_interval.data.name, gene_interval2.data.name, compute_discordance(gene_interval.data.sub_intervals, gene_interval2.data.sub_intervals))
        # for transcript_interval in gene_interval.data.sub_intervals:
        #     print(transcript_interval)
        #     for gene_interval2 in anno2_forward:
        #         for transcript_interval2 in gene_interval2:
        #             compute_discordance(transcript_interval, transcript_interval2)


