#!/usr/bin/env python3
import argparse
import pandas as pd
import os,sys
import time
import numpy as np
import re

parser = argparse.ArgumentParser(description="Drop features from fusion feature table")
parser.add_argument("-i", "--input", required=True, help="Input fusion feature table")
parser.add_argument("-o", "--output", required=True, help="Output fusion feature table with dropped features")
parser.add_argument("-d", "--drop_features", help="Comma-separated list of feature patterns to drop",
    default="burden")
parser.add_argument("-k", "--keep_features", help="Comma-separated list of feature patterns to keep")
args = parser.parse_args()
assert os.path.exists(args.input), f"Input file {args.input} does not exist."
assert not os.path.exists(args.output), f"Output file {args.output} already exists."

drop = args.drop_features.split(",")
print(f"Dropping features: {drop}")
# get keep features if specified
if args.keep_features:
    keep = args.keep_features.split(",")
    print(f"Keeping features: {keep}")
t = time.time()
df = pd.read_csv(args.input, sep="\t")
print(f"Read input file {args.input} in {time.time() - t:.2f} seconds.")
cols = df.columns.tolist()
print(f"Columns in input file: {cols}")
cols_keep = []
for c in cols:
    r=True
    # check if column matches any drop pattern
    if any(re.search(d, c) for d in drop):
        print(f"Dropping column: {c}")
        r = False
    # check if column matches any keep pattern
    # keep takes precedence over drop
    if args.keep_features:
        if any(re.search(k, c) for k in keep):
            print(f"Keeping column: {c}")
            r = True
    if r:
        cols_keep.append(c)
    
df = df[cols_keep]
print(f"Columns after filtering: {df.columns.tolist()}")

# Coerce counts to int
for c in df.columns:
    if 'count' in c:
        df[c] = df[c].astype(np.int64)  # Explicitly convert to integer type
t = time.time()
df.to_csv(args.output, sep="\t", index=False)
print(f"Wrote output file {args.output} in {time.time() - t:.2f} seconds.")
