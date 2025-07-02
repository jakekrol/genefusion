#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import argparse

# input: i) giggle file ii) column number of fileid iii) path to fid2type table

parser = argparse.ArgumentParser(description='Add specimen type to giggle file')
parser.add_argument('-i', '--input', type=str, required=True, help='Input giggle file')
parser.add_argument('-s', '--sample_column', type=int, required=True, help='Sample column index')
parser.add_argument('-l', '--lookup', type=str, default='/data/jake/genefusion/data/meta/fileid2sampletype.tsv', help='Path to file id to sample type mapping')
parser.add_argument('-o', '--output', type=str, required=True, help='Output giggle file with specimen type')
args = parser.parse_args()

args.sample_column -= 1

df_g = pd.read_csv(args.input, sep='\t', header=None)
df_l = pd.read_csv(args.lookup, sep='\t')
df_l.rename(columns={'File_ID': int(args.sample_column)}, inplace=True)
df_g = pd.merge(df_g, df_l, how='left', on=args.sample_column)
df_g.to_csv(args.output, sep='\t', header=None, index=False)
