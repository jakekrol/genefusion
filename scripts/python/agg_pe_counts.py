#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

# skiprows=1 to skip the header
SKIP=1
parser = argparse.ArgumentParser(description='Aggregate PE counts from multiple files')
parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input directory of gene-wise fusion count files')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for aggregated PE counts')
args = parser.parse_args()

assert not os.path.exists(args.output), f'Output file {args.output} already exists. Exiting.'   

files = os.listdir(args.input_dir)
files = [os.path.join(args.input_dir, f) for f in files]
n = len(files)
i = 1
# write all the files to the output file, exclude SKIP rows
t = time.time()
with open(args.output, 'a') as out:
    for f in files:
        with open(f, 'r') as fin:
            # skip the first SKIP lines
            for _ in range(SKIP):
                next(fin)
            # read the first line to get the header
            for line in fin:
                out.write(line)
        print(f'Processed {i}/{n} files')
        i += 1
print(f'Finished aggregating PE counts in {args.output} in {time.time()-t:.2f} seconds')

# now we add pe_count for gene duplicates
t = time.time()
df = pd.read_csv(args.output, sep='\t', header=None)
print(f'Loaded aggregated PE counts in {time.time()-t:.2f} seconds')
df.columns = ['left', 'right', 'pe_count']
t = time.time()
df = df.groupby(['left', 'right'])['pe_count'].sum().reset_index()
print(f'Combined PE counts for duplicate genes in {time.time()-t:.2f} seconds')
t = time.time()
df.to_csv(args.output, sep='\t', header=True, index=False)
print(f'Wrote aggregated PE counts to {args.output} in {time.time()-t:.2f} seconds')


    

            

