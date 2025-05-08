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

# goal: compute a null distribution of silhouette scores
# and compare it to the observed silhouette score
# of breakpoints
GENEFILE='/data/jake/genefusion/data/gene_file.txt.latest'
parser = argparse.ArgumentParser(description='Cluster breakpoints')
parser.add_argument('-i', '--input', type=str, required=True, help='Input file with breakpoints')
parser.add_argument('-l', '--left', help='left gene name', type=str, required=True)
parser.add_argument('-r', '--right',help='right gene name', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for clustered breakpoints')
parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed for reproducibility')
parser.add_argument('-p', '--permutations', type=int, default=100, help='Number of permutations for null distribution')
args = parser.parse_args()

# lookup start and end positions of genes
df_genes = pd.read_csv(GENEFILE, sep='\t', header=None)
df_genes.columns = ['chrm', 'start', 'end', 'name', 'strand']
l_start = df_genes[df_genes['name'] == args.left]['start'].values[0]
l_end = df_genes[df_genes['name'] == args.left]['end'].values[0]
r_start = df_genes[df_genes['name'] == args.right]['start'].values[0]
r_end = df_genes[df_genes['name'] == args.right]['end'].values[0]
# breakpoint data
df = pd.read_csv(args.input, sep='\t', header=None)
df.columns = [
        'hit_chrm', 'hit_start', 'hit_end', 'hit_name', 'hit_strand',
        'left_chrm', 'left_start', 'left_end', 'left_strand',
        'right_chrm', 'right_start', 'right_end', 'right_strand',
        'xxx', 'sample'
    ]
df = df[df['hit_name'] == args.right]
m = df.shape[0]
assert m > 0, f'No hits found for {args.right} in {args.input}'
df['left_mid'] = ((df['left_start'] + df['left_end']) / 2).apply(math.floor)
df['right_mid'] = ((df['right_start'] + df['right_end']) / 2).apply(math.floor)

# ignore large choice of k
K_CONSTRAINT = math.floor(m / 2)

# random data
random.seed(args.seed)

# get null distribution of silhouette scores
ss_null = []
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
    max_ss = -1
    for k in range(2, K_CONSTRAINT):
        kmeans = KMeans(n_clusters=k, random_state=args.seed)
        kmeans.fit(df_null)
        labels = kmeans.labels_
        ss = silhouette_score(df_null, labels)
        max_ss = max(max_ss, ss)
    ss_null.append(max_ss)

# compute observed silhouette score
max_ss = -1
for k in range(2, K_CONSTRAINT):
    kmeans = KMeans(n_clusters=k, random_state=args.seed)
    kmeans.fit(df[['left_mid', 'right_mid']])
    labels = kmeans.labels_
    ss = silhouette_score(df[['left_mid', 'right_mid']], labels)
    max_ss = max(max_ss, ss)
print(f'Observed silhouette score: {max_ss}')
    
# do histogram of null distribution and vline of observed silhouette score
plt.hist(ss_null, bins=20, alpha=0.5, color='blue', label='Null distribution')
plt.axvline(x=max_ss, color='red', linestyle='dashed', linewidth=2, label='Observed silhouette score')
plt.xlabel('Silhouette score')
plt.ylabel('Frequency')
plt.title('Null distribution of silhouette scores')
plt.legend()
plt.savefig(args.output)

    
    


