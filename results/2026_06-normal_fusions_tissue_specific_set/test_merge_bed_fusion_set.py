#!/usr/bin/env python
from polymerization.stix2fusion import *
from polymerization.io import *
import argparse

parser = argparse.ArgumentParser(description='Test merging of BED files for fusion sets.')
parser.add_argument('--bed', type=str, default='./grch37.genes.added.bed', help='Path to the input BED file')
parser.add_argument('--fusion_set', type=str, default='recurrent_normal_fusion_set.filtered.nohead.c1c2.tsv', help='Path to the input fusion set file')
args = parser.parse_args()

df_bed = read_bed(args.bed, gene_col_idx=3, uppercase=True)
df_fusion_set = read_fusion_set(args.fusion_set)
df_fusion_set = df_fusion_set.apply(lambda x: x.str.upper())
df_fusion_set_sort = left_sort_fusion_set(df_fusion_set, df_bed)
df_merge = merge_fusion_set_with_bed(df_fusion_set_sort, df_bed)

