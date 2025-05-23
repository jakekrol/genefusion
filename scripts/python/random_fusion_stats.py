#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

parser = argparse.ArgumentParser(description='Get stats for random fusions')
parser.add_argument('-q', '--query', type=str, required=True, help='Input query file (must be sorted s.t. left gene is col 1)')
parser.add_argument('-t', '--stats_tbl', type=str, required=True, help='Input fusion stats table')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for random fusion stats')
parser.add_argument('-n', '--n_samples', type=int, required=True, help='Number of random fusions to sample')
parser.add_argument('-s', '--allow_self', action='store_true', help='Allow self-fusions')   
parser.add_argument('-r', '--random_seed', type=int, default=42, help='Random seed for sampling')
args = parser.parse_args()

# load query
t = time.time()
df_q = pd.read_csv(args.query, sep='\t')
print(f'Loaded query file in {time.time() - t:.2f} seconds')
assert ('left' in df_q.columns) and ('right' in df_q.columns) and ('sort_result' in df_q.columns), f'Query file {args.query} should have "left", "right", and "sort_result" columns'

# load fusion stats
t= time.time()
df_f = pd.read_csv(args.stats_tbl, sep='\t')
print(f'Loaded fusion stats file in {time.time() - t:.2f} seconds')
assert ('left' in df_f.columns) and ('right' in df_f.columns), f'Fusion stats file {args.stats_tbl} should have "left" and "right" columns'

# drop the queries from the fusion stats
t= time.time()
df_f_filt = df_f.merge(df_q, on=['left', 'right'], how='left', indicator=True)
df_f_filt = df_f_filt[df_f_filt['_merge'] == 'left_only'].drop(columns=['_merge'])
print(f'Dropped queries from fusion stats in {time.time() - t:.2f} seconds')

# remove self-fusions if not allowed
t= time.time()
if not args.allow_self:
    df_f_filt = df_f_filt[df_f_filt['left'] != df_f_filt['right']]
print(f'Removed self-fusions in {time.time() - t:.2f} seconds')

# sample random fusions
t = time.time()
df_random = df_f_filt.sample(n=args.n_samples, random_state=args.random_seed)
df_random.to_csv(args.output, sep='\t', index=False)


