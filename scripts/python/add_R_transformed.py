#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time
import numpy as np

parser = argparse.ArgumentParser(description='Add |R-1| column to gene fusion stats file')
parser.add_argument('-i', '--input', type=str, required=True, help='Input gene fusion stats file with R column')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file with R transformed column')
parser.add_argument('-t', '--type', default='tumor', choices=['tumor', 'normal'], help='Type of fusion table (tumor or normal)')
args = parser.parse_args()

if os.path.exists(args.output):
    print(f'Output file {args.output} already exists. Exiting.')
    sys.exit(1)

t = time.time()
df = pd.read_csv(args.input, sep='\t')
print(f'Input read in {time.time() - t:.2f} seconds')

# apply the transformation
t = time.time()
df[f'R_transformed_{args.type}'] = np.abs(df[f'R_{args.type}']-1) # nans will be preserved
print(f'R transformation applied in {time.time() - t:.2f} seconds')
t = time.time()
df.to_csv(args.output, sep='\t', index=False)
print(f'Output written in {time.time() - t:.2f} seconds')
