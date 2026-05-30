#!/usr/bin/env python3

import argparse
from polymerization.stix2fusion import *
from polymerization.io import *

parser = argparse.ArgumentParser(description='Sort fusion set by genomic position')
parser.add_argument('--fusion_set', default='./query_fusions.tsv', help='path to fusion set file')
parser.add_argument('--bedfile', default='./grch37.genes.bed.added', help='path to bed file')
parser.add_argument('--outfile', default='./query_fusions_sorted.tsv', help='path to output sorted fusion set file')
args = parser.parse_args()

df_fusion = read_fusion_set(args.fusion_set, reader='pandas')
df_bed = read_bed(args.bedfile, gene_col_idx=3)
df_fusion_sort = left_sort_fusion_set(df_fusion, df_bed)
df_fusion_sort.to_csv(args.outfile, sep='\t', index=False)