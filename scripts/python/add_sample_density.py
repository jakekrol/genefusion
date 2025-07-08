#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

parser = argparse.ArgumentParser(description='Compute and add sample density to fusion table')
parser.add_argument('-i', '--input', type=str, required=True, help='Input fusion table file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for fusion table with sample density')
parser.add_argument('-t', '--type', default='tumor', choices=['tumor', 'normal'], help='Type of fusion table (tumor or normal)')
args = parser.parse_args()

t = time.time()
df = pd.read_csv(args.input, sep='\t')
df[f'pe_count_{args.type}'] = df[f'pe_count_{args.type}'].astype(int)
df[f'sample_count_{args.type}'] = df[f'sample_count_{args.type}'].astype(int)
print(f'Loaded fusion table file in {time.time() - t:.2f} seconds')
t = time.time()
df[f'sample_density_{args.type}'] = df[f'sample_count_{args.type}'] / df[f'pe_count_{args.type}']
print(f'Computed sample density in {time.time() - t:.2f} seconds')
t = time.time()
df.to_csv(args.output, sep='\t', index=False)
print(f'Wrote result to {args.output} in {time.time() - t:.2f} seconds')
