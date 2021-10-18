#!/usr/bin/env python

from __future__ import print_function, division
import argparse

from intervaltree import Interval, IntervalTree

from classify_annotation import gff_to_intervals

def compute_overlap(anno1: IntervalTree, anno2: IntervalTree) -> float:
    overall_length = 0
    overall_overlap = 0
    for interval in anno1:
        # note: intervals are "half open", i.e. begin is in the interval, end is outside
        # by contrast GFF3 intervals are open
        interval_length = interval.length()
        overall_length += interval_length
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--use_exon', default=False, action='store_true')
    parser.add_argument('anno1_file', type=argparse.FileType())
    parser.add_argument('anno2_file', type=argparse.FileType())
    args = parser.parse_args()

    (anno1_forward, anno1_reverse) = gff_to_intervals(args.anno1_file)
    (anno2_forward, anno2_reverse) = gff_to_intervals(args.anno2_file)
    for gene_interval in anno1_forward:
        print(gene_interval)
        break


