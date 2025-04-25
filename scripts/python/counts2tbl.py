#!/usr/bin/env python3

# input: a line delimited list of paths to a sample's fusion count files
# for ex just pipe the grep'd sample name 
# output: genefusion count table for  a sample


import argparse
import os, sys
import pandas as pd

# expect the fileformat of sample wise counts to be:
# sample.gene.(v ...).chromosome.strandd.left.right
# where (v ...) are optional attributes such as version num.
N_EXPECTED_FILE_ATTR = 6
parser = argparse.ArgumentParser(description='Convert gene fusion count files to a table.')
parser.add_argument('-i', '--input', type=str, help='Input directory')
parser.add_argument('-i', '--output', type=str, help='Output file')
args = parser.parse_args()

# get files
files = os.listdir(args.input)
files = [os.path.join(args.input, f) for f in files if os.path.isfile(os.path.join(args.input, f))]
assert len(files) > 0, "No files provided. Please provide a list of gene fusion count files as stdin."

l = []
n = len(files)
i = 0
# build  list of 3-tuples (gene_i, gene_j, count)
for file in files:
    # careful genes often have '.' in their names
    fname = os.path.basename(file)
    g_i = fname.split('.')[1]
    d = len(fname.split('.')) - N_EXPECTED_FILE_ATTR
    if d > 0:
        for i in range(2,2+d):
            g_i += '.' + fname.split('.')[i]
    with open(file) as f:
        for line in f:
            g_j,count = f.readline().strip().split('\t')
            l.append((g_i, g_j, count))
    print(f'Processed {i+1}/{n} files')
df = pd.DataFrame(l, columns=['gene_i', 'gene_j', 'count'])
df['count'] = df['count'].astype(int)
df.to_csv(args.output, sep='\t', index=False)
