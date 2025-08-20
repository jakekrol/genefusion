#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

parser = argparse.ArgumentParser(description='Join stats onto queries')
parser.add_argument('-q', '--query', type=str, required=True, help='Input query file (must be sorted s.t. left gene is col 1)')
parser.add_argument('-t', '--tbl', type=str, required=True, help='Input fusion table file')
parser.add_argument('-b', '--burden', type=str, required=True, help='Input burden table file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file')  
parser.add_argument('-d', '--drop', action='store_true', help='Drop fusions with sort_result=-1')   
args = parser.parse_args()

# load query
t = time.time()
df_q = pd.read_csv(args.query, sep='\t')
print(f'Loaded query file in {time.time() - t:.2f} seconds')
assert ('left' in df_q.columns) and ('right' in df_q.columns), f'Query file {args.query} should have "left" and "right" columns'

# load fusion
t= time.time()
df_f = pd.read_csv(args.tbl, sep='\t')
assert ('left' in df_f.columns) and ('right' in df_f.columns), f'Fusion table {args.tbl} should have "left" and "right" columns'    
print(f'Loaded fusion table file in {time.time() - t:.2f} seconds')

# load burden
t= time.time()
df_b = pd.read_csv(args.burden, sep='\t',header=None)
df_b.columns = ['gene', 'burden_total']
print(f'Loaded burden table file in {time.time() - t:.2f} seconds')

# merge
t= time.time()
df_merged = pd.merge(df_q, df_f, on=['left', 'right'], how='left')
print(f'Merged query and fusion table in {time.time() - t:.2f} seconds')

# write
t= time.time()
df_merged.to_csv(args.output, sep='\t', index=False)
print(f'Wrote output file in {time.time() - t:.2f} seconds')

# handle case where sort_result is 1, yet no PE evidence in the fusion table
# if sort_result is 1 and NaN, then fill pe_count and sample_count with 0
t = time.time()
mask = (df_merged['sort_result'] == 1) & (df_merged['pe_count'].isna())
df_merged.loc[mask, 'pe_count'] = 0
df_merged.loc[mask, 'sample_count'] = 0
print(f'Filled successful sort_result NaN in {time.time() - t:.2f} seconds')
# add burden to these rows
t = time.time()
for i,val in mask.items():
    if val:
        left = df_merged.loc[i, 'left']
        right = df_merged.loc[i, 'right']
        burden_l = df_b[df_b['gene'] == left]['burden_total'].values[0]
        burden_r = df_b[df_b['gene'] == right]['burden_total'].values[0]
        df_merged.loc[i, 'burden_total_left'] = burden_l
        df_merged.loc[i, 'burden_total_right'] = burden_r
print(f'Added burden in {time.time() - t:.2f} seconds')
# drop fusions with sort_result=-1
if args.drop:
    t= time.time()
    df_merged = df_merged[df_merged['sort_result'] != -1]
    print(f'Dropped fusions with sort_result=-1 in {time.time() - t:.2f} seconds')
# write again
t= time.time()
df_merged.to_csv(args.output, sep='\t', index=False)
print(f'Wrote output file with sort_result=1 NAs filled in {time.time() - t:.2f} seconds')    

