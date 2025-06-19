#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import math
import random
from genefusion.genefusion import clark_evans_R
import seaborn as sns

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
parser.add_argument('-l', '--left', help='left gene name', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True, help='Output table path')
parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed for reproducibility')
parser.add_argument('-p', '--permutations', type=int, default=1000, help='Number of permutations for null distribution')
parser.add_argument('--min_n', type=int, default=10, help='Minimum number of breakpoints for a fusion to be considered')
args = parser.parse_args()
random.seed(args.seed)

# lookup start and end positions of genes
df_genes = pd.read_csv(GENEFILE, sep='\t', header=None)
df_genes.columns = ['chrm', 'start', 'end', 'name', 'strand']
l_start = df_genes[df_genes['name'] == args.left]['start'].values[0]
l_end = df_genes[df_genes['name'] == args.left]['end'].values[0]
# breakpoint data
df = pd.read_csv(args.input, sep='\t', header=None)
df.columns = [
        'hit_chrm', 'hit_start', 'hit_end', 'hit_name', 'hit_strand',
        'left_chrm', 'left_start', 'left_end', 'left_strand',
        'right_chrm', 'right_start', 'right_end', 'right_strand',
        'xxx', 'sample', 'specimen'
    ]
output = {'right': [], 'R': [], 'p_value': []}
n_rights = df['hit_name'].nunique()
counter = 1
for value, subset in df.groupby('hit_name'):
    m = subset.shape[0]
    if m < args.min_n:
        print(f'Skipping {value} with only {m} breakpoints')
        counter += 1
        continue
    r_start = df_genes[df_genes['name'] == value]['start'].values[0]
    r_end = df_genes[df_genes['name'] == value]['end'].values[0]
    subset['left_mid'] = ((subset['left_start'] + subset['left_end']) / 2).apply(math.floor)
    subset['right_mid'] = ((subset['right_start'] + subset['right_end']) / 2).apply(math.floor)

    # compute R stat for breakpoint randomness
    x = subset['left_mid'].values
    y = subset['right_mid'].values
    # emprical R stat
    R_emp = clark_evans_R(x, y, x_start=l_start, x_end=l_end, y_start=r_start, y_end=r_end)['R']

    # get null distribution of R values
    R_null = []
    for i in range(args.permutations):
        x = []
        a = l_start
        b = l_end
        for i in range(m):
            x.append(random.randint(a, b))
        y = []
        a = r_start
        b = r_end
        for i in range(m):
            y.append(random.randint(a, b))
        df_null = pd.DataFrame({'x': x, 'y': y})
        # get R stat for null distribution
        rand = clark_evans_R(df_null['x'].values, df_null['y'].values, x_start=l_start, x_end=l_end, y_start=r_start, y_end=r_end)
        R_null.append(rand['R'])

    # compute the monte carlo p-value by two-tailed test
    R_null = np.array(R_null)
    mean_null = np.mean(R_null)
    # test whether a null value is as or more extreme than the observed value
    def is_null_extreme(null, observed, mean_null):
        return abs(null - mean_null) >= abs(observed - mean_null)
    numerator = 0
    for x_i in R_null:
        if is_null_extreme(x_i, R_emp, mean_null):
            numerator += 1
    p_value = (numerator + 1) / (args.permutations + 1)
    output['right'].append(value)
    output['R'].append(R_emp)
    output['p_value'].append(p_value)
    print(f'Processed {counter}/{n_rights} fusions: {value} | R: {R_emp:.3f} | p-value: {p_value:.3f}')
    counter += 1
output_df = pd.DataFrame(output)
output_df['left'] = args.left
output_df.to_csv(args.output, sep='\t', index=False)
print(f'Output written to {args.output}')
