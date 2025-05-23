#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

parser = argparse.ArgumentParser(description='Compute and add sample density to fusion table')
parser.add_argument('-i', '--input', type=str, required=True, help='Input fusion table file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for fusion table with sample density')
args = parser.parse_args()

t = time.time()
df = pd.read_csv(args.input, sep='\t')
df['pe_count'] = df['pe_count'].astype(int)
df['sample_count'] = df['sample_count'].astype(int)
print(f'Loaded fusion table file in {time.time() - t:.2f} seconds')
t = time.time()
df['sample_density'] = df['sample_count'] / df['pe_count']
print(f'Computed sample density in {time.time() - t:.2f} seconds')
t = time.time()
df.to_csv(args.output, sep='\t', index=False)
print(f'Wrote result to {args.output} in {time.time() - t:.2f} seconds')
