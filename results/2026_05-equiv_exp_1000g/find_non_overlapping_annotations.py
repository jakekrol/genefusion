#!/usr/bin/env python3
import argparse
import pandas as pd
from polymerization.io import read_fusion_set, read_bed
from polymerization.stix2fusion import merge_fusion_set_with_bed, left_sort_fusion_set

parser=argparse.ArgumentParser(description='Find gene pairs with non-overlapping annotations')
parser.add_argument('--fusion_set', default='highest_read_depth_fusion_genepairs.tsv')
parser.add_argument('--gene_bedfile', default='../2025_04-gene_bedfile_cln/grch37.genes.bed')
parser.add_argument('--out_file', default='fusion_genepairs_non_overlapping_annotations.tsv')
args=parser.parse_args()

df_fusion = read_fusion_set(args.fusion_set, sep='\t')
df_bed = read_bed(args.gene_bedfile, sep='\t', gene_col_idx=3)
df_fusion_sort = left_sort_fusion_set(df_fusion, df_bed)
df_merge = merge_fusion_set_with_bed(df_fusion_sort, df_bed, merger='pandas')
def check_overlap(chr_left, start_left, end_left, chr_right, start_right, end_right):
    if chr_left != chr_right:
        return False
    return not (end_left <= start_right or end_right <= start_left)
for idx, row in df_merge.iterrows():
    overlap = check_overlap(row['chromosome_left'], row['start_left'], row['end_left'], row['chromosome_right'], row['start_right'], row['end_right'])
    df_merge.at[idx, 'overlap'] = overlap
df_non_overlapping = df_merge[df_merge['overlap'] == False]
df_non_overlapping.to_csv(args.out_file, sep='\t', index=False)


