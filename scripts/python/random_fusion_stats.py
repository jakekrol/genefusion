#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time
import subprocess

parser = argparse.ArgumentParser(description='Get stats for random fusions')
parser.add_argument('-q', '--query', type=str, required=True, help='Input query file (must be sorted s.t. left gene is col 1)')
parser.add_argument('-t', '--stats_tbl', type=str, required=True, help='Input fusion stats table')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for random fusion stats')
parser.add_argument('-n', '--n_samples', type=int, help='Number of random fusions to sample')
parser.add_argument('-s', '--allow_self', action='store_true', help='Allow self-fusions')   
parser.add_argument('-r', '--random_seed', type=int, default=42, help='Random seed for sampling')
args = parser.parse_args()

    
# load query
t = time.time()
df_q = pd.read_csv(args.query, sep='\t')
print(f'Loaded query file in {time.time() - t:.2f} seconds')
assert ('left' in df_q.columns) and ('right' in df_q.columns), f'Query file {args.query} should have "left" and "right" columns'

# load fusion stats
t= time.time()
df_f = pd.read_csv(args.stats_tbl, sep='\t')
print(f'Loaded fusion stats file in {time.time() - t:.2f} seconds')
assert ('left' in df_f.columns) and ('right' in df_f.columns), f'Fusion stats file {args.stats_tbl} should have "left" and "right" columns'

# drop the queries from the fusion stats
t= time.time()
df_f_filt = df_f.merge(df_q, on=['left', 'right'], how='left', indicator=True) # add indicator column to see which rows are from the query
df_f_filt = df_f_filt[df_f_filt['_merge'] == 'left_only'].drop(columns=['_merge']) # drop the rows that are in the query
# keep only the columns from the fusion stats
cols_to_keep = [col for col in df_f_filt.columns if col.endswith('_x') or col in ['left', 'right']]
df_f_filt = df_f_filt[cols_to_keep]
df_f_filt.columns = [col[:-2] if col.endswith('_x') else col for col in df_f_filt.columns]
print(f'Dropped queries from fusion stats in {time.time() - t:.2f} seconds')

# remove self-fusions if not allowed
t= time.time()
if not args.allow_self:
    df_f_filt = df_f_filt[df_f_filt['left'] != df_f_filt['right']]
print(f'Removed self-fusions in {time.time() - t:.2f} seconds')

# drop duplicates and permutations
# write tmp
tmp_file = args.stats_tbl + '.tmp'
t = time.time()
df_f_filt.to_csv(tmp_file, sep='\t', index=False)
print(f'Wrote temp file {tmp_file} in {time.time() - t:.2f} seconds')
ncols = len(df_f_filt.columns)
s = ''
for i in range(3,ncols+1):
    s += f'${i}, '
s = s[:-2]  # Remove trailing comma and space
cmd_str = f"awk 'BEGIN {{OFS=\"\\t\"}} {{print $1, $2, $1 \"zzz\" $2, $2 \"zzz\" $1, {s}}}' {tmp_file} > {args.stats_tbl}.tmp.keyed"
print("Running command: " + cmd_str)
subprocess.run(cmd_str, shell=True, check=True)
print(f'Created keyed fusion stats file in {time.time() - t:.2f} seconds')
t=time.time()
df_keyed = pd.read_csv(args.stats_tbl + '.tmp' + '.keyed', sep='\t')
print(f'Loaded keyed fusion stats file in {time.time() - t:.2f} seconds')
# drop duplicates
t= time.time()
df_keyed = df_keyed.drop_duplicates(subset=['left', 'right'])
df_keyed = df_keyed.drop(columns=['leftzzzright', 'rightzzzleft'])
print(f'Dropped duplicates in {time.time() - t:.2f} seconds')

df_f_filt = df_keyed

# fill nas
cols = df_f_filt.columns.tolist()
fill0 = [c for c in cols if ('count' in c) or ('transformed' in c)]
for c in fill0:
    df_f_filt[c] = df_f_filt[c].fillna(0)
fill1 = [c for c in cols if 'density' in c]
for c in fill1:
    df_f_filt[c] = df_f_filt[c].fillna(1)

# sample random fusions
t = time.time()
if not args.n_samples:
    print(f'No n_samples specified, will include all {df_f_filt.shape[0]} non-query fusions in output')
    df_random = df_f_filt
    df_random.to_csv(args.output, sep='\t', index=False)
    print(f'Wrote all non-query fusions in {time.time() - t:.2f} seconds')
else:
    df_random = df_f_filt.sample(n=args.n_samples, random_state=args.random_seed)
    df_random.to_csv(args.output, sep='\t', index=False)
    print(f'Sampled and wrote {args.n_samples} random fusions in {time.time() - t:.2f} seconds')


