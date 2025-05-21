#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import re
import time

# input: giggle file w/ sample column
# output: giggle file w/ cleaned sample column
# note: for intersect files, sample column is 15 (1-indexed)

parser = argparse.ArgumentParser(description="Clean sample names in giggle file")
parser.add_argument('-i', '--input', type=str, required=True, help='Input giggle file')
parser.add_argument('-s', '--sample_col', type=str, required=True, help='Sample column index')
parser.add_argument('-o', '--output', type=str, required=True, help='Output giggle file')
args = parser.parse_args()

args.sample_col = int(args.sample_col) - 1  # convert to 0-indexed

def cln_sample_names(x):
    # remove index prefix
    x = os.path.basename(x)
    # keep only the fileid
    x = x.split(".")[0]
    return x

df = pd.read_csv(args.input, sep="\t", header=None)
t = time.time()
df[args.sample_col] = df[args.sample_col].apply(cln_sample_names)
df.to_csv(args.output, sep="\t", header=None, index=False)
print(f'Cleaned sample col of {args.input} in {time.time() - t:.2f} seconds to {args.output}')
