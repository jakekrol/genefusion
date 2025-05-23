#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

parser = argparse.ArgumentParser(description='Count fusions in query file')
parser.add_argument('-i', '--input', type=str, required=True, help='Input giggle file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
parser.add_argument('-l', '--left', type=str, help='Left gene name')
parser.add_argument('-z', '--hack', action='store_true', help='Hack to parse left gene name from the input file name')
parser.add_argument('-r', '--right_gene_col_idx', type=int, default=4, help='Right gene column index (1-indexed)')

args = parser.parse_args()
if args.hack:
    args.left = '.'.join(os.path.basename(args.input).split('.')[:-4]) # remove the last 4 elements (chrom, strand, start, end)
    assert len(args.left) > 0, f'Left gene name could not be parsed from the input file name: result={args.left}. Please provide a left gene name using -l option.'
    print(f'Parsed left gene name from input file name: {args.left}')

args.right_gene_col_idx -= 1

df = pd.read_csv(args.input, sep='\t', header=None, usecols=[args.right_gene_col_idx])
df.columns = ['right']
s = df['right'].value_counts()
df = s.to_frame()
df.reset_index(inplace=True) # make right a column
df['left'] = args.left
df = df.rename(columns={'count': 'pe_count'})
df = df[['left', 'right', 'pe_count']]
df.sort_values(by=['left', 'right'], inplace=True)
df.to_csv(args.output, sep='\t', header=True,index=False)
