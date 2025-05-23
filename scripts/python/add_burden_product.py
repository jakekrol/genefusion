#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time
import numpy as np

parser = argparse.ArgumentParser(description='Compute and add burden product to fusion table')
parser.add_argument('-i', '--input', type=str, required=True, help='Input fusion table file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for fusion table with burden product')
args = parser.parse_args()

t = time.time()
df = pd.read_csv(args.input, sep='\t')
print(f'Loaded fusion table file in {time.time() - t:.2f} seconds')
df['burden_total_left'] = df['burden_total_left'].astype(int) 
df['burden_total_right'] = df['burden_total_right'].astype(int)
t = time.time()
df['burden_product'] = df['burden_total_left'] * df['burden_total_right']
print(f'Computed burden product in {time.time() - t:.2f} seconds')
t= time.time()
df['log10_burden_product'] = df['burden_product'].apply(lambda x: 0 if x == 0 else round(np.log10(x), 3))
print(f'Computed burden product in {time.time() - t:.2f} seconds')
t = time.time()
df.to_csv(args.output, sep='\t', index=False)
print(f'Wrote result to {args.output} in {time.time() - t:.2f} seconds')
