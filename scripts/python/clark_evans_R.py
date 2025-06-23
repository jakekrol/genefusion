#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import math
import numpy as np
from genefusion.genefusion import clark_evans_R
import time

# goal: compute R stat and p-value for all fusions in a giggle file
# key inputs
# input: i) an intersected giggle file with left and right breakpoint data of fusions
# input: ii) left gene name
# key outputs
# output: i) R value
# output: ii) p-value
GENEFILE='/data/jake/genefusion/data/gene_file.txt.latest'
parser = argparse.ArgumentParser(description='Cluster breakpoints')
parser.add_argument('-i', '--input', type=str, required=True, help='Input file with breakpoints')
parser.add_argument('-l', '--left', help='left gene name', type=str)
parser.add_argument('-o', '--output', type=str, required=True, help='Output table path')
# parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed for reproducibility')
# parser.add_argument('-p', '--permutations', type=int, default=1000, help='Number of permutations for null distribution')
parser.add_argument('-n', '--min_n', type=int, default=10, help='Minimum number of breakpoints for a fusion to be considered')
parser.add_argument('-z', '--hack', action='store_true', help='Hack to parse left gene name from the input file name')
parser.add_argument('-d', '--downsample', type = int, default=-1, help='Downsample the targed grouped df to this number of rows. If -1, no downsampling is performed.')
parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility (default: 42)')
args = parser.parse_args()
args = parser.parse_args()
if args.hack:
    args.left = '.'.join(os.path.basename(args.input).split('.')[:-4]) # remove the last 4 elements (chrom, strand, start, end)
    assert len(args.left) > 0, f'Left gene name could not be parsed from the input file name: result={args.left}. Please provide a left gene name using -l option.'
    print(f'Parsed left gene name from input file name: {args.left}')
print(args)

# random.seed(args.seed)

# lookup start and end positions of genes
df_genes = pd.read_csv(GENEFILE, sep='\t', header=None)
df_genes.columns = ['chrm', 'start', 'end', 'name', 'strand']
l_start = df_genes[df_genes['name'] == args.left]['start'].values[0]
l_end = df_genes[df_genes['name'] == args.left]['end'].values[0]
# breakpoint data
t = time.time()
# - necessary columns: (1, hit_start), (2, hit_end), (3, hit_name), (5, left_start) (6, left_end), (10, right_start), (11, right_end)
# - excluded columns: (0, hit_chrm), (4, hit_strand), (5, left_chrm), (8, left_strand), (9, right_chrm), (12, right_strand), (13, xxx) (14, sample) (15, specimen)
df = pd.read_csv(args.input, sep='\t', header=None, usecols=[1, 2, 3, 5, 6, 10, 11])
print(f'Loaded input file in {time.time() - t:.2f} seconds')
df.columns = [
        'hit_start', 'hit_end', 'hit_name',
        'left_start', 'left_end',
        'right_start', 'right_end'
    ]
# output = {'right': [], 'R': [], 'p_value': []}
output = {'right': [], 'R': []}
n_rights = df['hit_name'].nunique()
t = time.time()
for value, subset in df.groupby('hit_name'):
    m = subset.shape[0]
    if (args.downsample > 0) and (m > args.downsample):
        subset = subset.sample(n=args.downsample, random_state=args.seed)
    m = subset.shape[0]
    if m < args.min_n:
        continue
    r_start = df_genes[df_genes['name'] == value]['start'].values[0]
    r_end = df_genes[df_genes['name'] == value]['end'].values[0]
    subset['left_mid'] = np.floor( (subset['left_start'] + subset['left_end']) / 2).astype(int)
    subset['right_mid'] = np.floor( (subset['right_start'] + subset['right_end']) / 2).astype(int)

    # compute R stat for breakpoint randomness
    x = subset['left_mid'].values
    y = subset['right_mid'].values
    # emprical R stat
    R_emp = clark_evans_R(x, y, x_start=l_start, x_end=l_end, y_start=r_start, y_end=r_end)['R']
    output['right'].append(value)
    output['R'].append(R_emp)
print(f'Processed all fusions in {time.time() - t:.2f} seconds')

    # monte carlo was too slow, so we are not using it

    # # get null distribution of R values
    # R_null = []
    # for i in range(args.permutations):
    #     x = []
    #     a = l_start
    #     b = l_end
    #     for i in range(m):
    #         x.append(random.randint(a, b))
    #     y = []
    #     a = r_start
    #     b = r_end
    #     for i in range(m):
    #         y.append(random.randint(a, b))
    #     df_null = pd.DataFrame({'x': x, 'y': y})
    #     # get R stat for null distribution
    #     rand = clark_evans_R(df_null['x'].values, df_null['y'].values, x_start=l_start, x_end=l_end, y_start=r_start, y_end=r_end)
    #     R_null.append(rand['R'])

    # # compute the monte carlo p-value by two-tailed test
    # R_null = np.array(R_null)
    # mean_null = np.mean(R_null)
    # # test whether a null value is as or more extreme than the observed value
    # def is_null_extreme(null, observed, mean_null):
    #     return abs(null - mean_null) >= abs(observed - mean_null)
    # numerator = 0
    # for x_i in R_null:
    #     if is_null_extreme(x_i, R_emp, mean_null):
    #         numerator += 1
    # p_value = (numerator + 1) / (args.permutations + 1)
    # output['right'].append(value)
    # output['R'].append(R_emp)
    # output['p_value'].append(p_value)
output_df = pd.DataFrame(output)
output_df['left'] = args.left
output_df = output_df[['left', 'right', 'R']]
t = time.time()
output_df.to_csv(args.output, sep='\t', index=False)
print(f'Saved output file to {args.output} in {time.time() - t:.2f} seconds')
