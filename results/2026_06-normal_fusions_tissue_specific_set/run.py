#!/usr/bin/env python3

import pandas as pd
import os,sys
import argparse

parser = argparse.ArgumentParser(description="make fusion set")
parser.add_argument('--input', default = '../../data/2026_06-normal_tissue_specific_fusions/recurrent_normal_tissue_specific.tsv', help='input file')
parser.add_argument('--output', default = 'recurrent_normal_fusion_set.tsv', help='output file')
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t')
col_keep = ['up_gene', 'dw_gene', 'sample source', 'rt_pcr_confirmed', 'western_blot_or_nmr_confirmed']
df = df[col_keep]
df['rt_pcr_confirmed'] = df['rt_pcr_confirmed'].fillna(0).astype(int)
df['western_blot_or_nmr_confirmed'] = df['western_blot_or_nmr_confirmed'].fillna(0).astype(int)
df.rename(columns={'sample source': 'sample_source'}, inplace=True)
df = df.drop_duplicates().reset_index(drop=True)
df.to_csv(args.output, sep='\t', index=False)
