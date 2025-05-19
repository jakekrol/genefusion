#!/usr/bin/env python3

import pandas as pd
import os
import sys
import argparse
import time

# ex: ./add_burden_col.py -f tumor_pe_and_sample_count.tsv -b burden_total_tumor.tsv -o tumor_pe_sample_and_burden.tsv -k1 0 -k2 0 -n burden_total_left
parser = argparse.ArgumentParser(description='Add burden column to gene fusion table')
parser.add_argument('-f', '--fusion_file', type=str, required=True, help='Input gene fusion table file')
parser.add_argument('-b', '--burden_file', type=str, required=True, help='Input burden file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for gene fusion table with burden')
parser.add_argument('-h1', '--fusion_header', action='store_true', help='Fusion file has header')
parser.add_argument('-h2', '--burden_header', action='store_true', help='Burden file has header')
parser.add_argument('-k1', '--fusion_key', type=int, required=True, help='Column index for fusion key in fusion file')
parser.add_argument('-k2', '--burden_key', type=int, required=True, help='Column index for burden key in burden file')
parser.add_argument('-n', '--burden_name', type=str, required=True, help='Column name for burden in output file')

args = parser.parse_args()
if os.path.exists(args.output):
    print(f'Output file {args.output} already exists. Exiting.')
    sys.exit(1)

# Read the fusion file
t = time.time()
if args.fusion_header:
    df_fusion = pd.read_csv(args.fusion_file, sep='\t')
else: 
    df_fusion = pd.read_csv(args.fusion_file, sep='\t', header=None)
print(f'Fusion file read in {time.time() - t:.2f} seconds')

# Read the burden file
t = time.time()
if args.burden_header:
    df_burden = pd.read_csv(args.burden_file, sep='\t')
else:
    df_burden = pd.read_csv(args.burden_file, sep='\t', header=None)
print(f'Burden file read in {time.time() - t:.2f} seconds')

# Merge the two dataframes on the fusion key
t = time.time()
# update the key if header is present
if args.fusion_header:
    args.fusion_key = df_fusion.columns[args.fusion_key]
if args.burden_header:
    args.burden_key = df_burden.columns[args.burden_key]
# print keys
print(f'Merging on fusion key: {args.fusion_key} and burden key: {args.burden_key}')
df_merged = pd.merge(df_fusion, df_burden, left_on=args.fusion_key, right_on=args.burden_key, how='left')
# remove the burden gene key from the merged dataframe
if args.burden_key in df_merged.columns:
    df_merged.drop(columns=[args.burden_key], inplace=True)
print(df_merged.head())
print(f'Merged file created in {time.time() - t:.2f} seconds')

# Write the merged file
print(f'df_fusion columns: {df_fusion.columns.tolist()}')
print('df_merged columns: ', df_merged.columns.tolist())
columns = df_fusion.columns.tolist() + [args.burden_name]
print(f'New columns: {columns}')
df_merged.columns = columns
t = time.time()
df_merged.to_csv(args.output, sep='\t', index=False)
print(f'Merged file saved in {time.time() - t:.2f} seconds')

