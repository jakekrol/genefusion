#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

parser = argparse.ArgumentParser(description='Randomly sample a two column fusion table')
parser.add_argument('-b', '--bedfile', type=str, default="/data/jake/genefusion/data/gene_file.txt.latest", help='Input gene bed file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
parser.add_argument('-n', '--num_samples', type=int, required=True, help='Number of samples to generate')
args = parser.parse_args()

df_b = pd.read_csv(args.bedfile, sep='\t', header=None, usecols=[3])
df_b.columns = ['gene']
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