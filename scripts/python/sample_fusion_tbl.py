#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

print("Deprecated: This script is deprecated.")
sys.exit(1)
parser = argparse.ArgumentParser(description='Randomly sample a two column fusion table')
parser.add_argument('-b', '--bedfile', type=str, default="/data/jake/genefusion/data/gene_file.txt.latest", help='Input gene bed file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
parser.add_argument('-n', '--num_samples', type=int, required=True, help='Number of samples to generate')
parser.add_argument('-x', '--exclude', type=str, default=None, help='Exclude fusions from sampling')
parser.add_argument('-hx', '--exclude_header', action='store_true', help='Exclude file has header')
args = parser.parse_args()

df_b = pd.read_csv(args.bedfile, sep='\t', header=None, usecols=[3])
df_b.columns = ['gene']
if args.exclude:
    if args.exclude_header:
        df_exclude = pd.read_csv(args.exclude, sep='\t')
    else:
        df_exclude = pd.read_csv(args.exclude, sep='\t', header=None)
for i in range(args.num_samples):
    # sample two genes
    left = df_b.sample(1).values[0][0]
    right = df_b.sample(1).values[0][0]
    # make sure they are not the same
    while left == right:
        right = df_b.sample(1).values[0][0]
    # write to file
    with open(args.output, 'a') as out:
        out.write(f'{left}\t{right}\n')
    print(f'Processed {i+1}/{args.num_samples} rows')
print(f'Wrote {args.num_samples} rows to {args.output}')