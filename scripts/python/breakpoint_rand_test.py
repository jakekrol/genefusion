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

# key inputs
# input: i) an intersected giggle file with left and right breakpoint data of fusions
# input: ii) left gene name
# input: iii) right gene name
# goal: compute R stat for breakpoint randomness
# optionally do a Monte Carlo test to get a p-value for the R stat
# key outputs
# output: i) a scatter plot of the breakpoints with probability distribution for each axis
# output: ii) a histogram of null distribution of R stats with vline of empirical R and p-value
GENEFILE='/data/jake/genefusion/data/gene_file.txt.latest'
parser = argparse.ArgumentParser(description='Cluster breakpoints')
parser.add_argument('-i', '--input', type=str, required=True, help='Input file with breakpoints')
parser.add_argument('-l', '--left', help='left gene name', type=str, required=True)
parser.add_argument('-r', '--right',help='right gene name', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True, help='Output file pattern')
parser.add_argument('--log', action='store_true', help='log10 transform the observed and expected distances returned by clark_evans_R')
parser.add_argument('-f', '--fontsize', type=int, default=18, help='Font size for the plot')
parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed for reproducibility')
parser.add_argument('-p', '--permutations', type=int, default=1000, help='Number of permutations for null distribution')
parser.add_argument('--plot', action='store_true', help='Plot the breakpoints with probability distribution for each axis')
parser.add_argument('--mc', action='store_true', help='Use Monte Carlo method for null distribution')
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
        'xxx', 'sample', 'specimen'
    ]
df = df[df['hit_name'] == args.right]
m = df.shape[0]
assert m > 0, f'No hits found for {args.right} in {args.input}'
df['left_mid'] = ((df['left_start'] + df['left_end']) / 2).apply(math.floor)
df['right_mid'] = ((df['right_start'] + df['right_end']) / 2).apply(math.floor)

# compute R stat for breakpoint randomness
x = df['left_mid'].values
y = df['right_mid'].values
# emprical R stat
emp = clark_evans_R(x, y, x_start=l_start, x_end=l_end, y_start=r_start, y_end=r_end, log10=args.log)
emp = {key: round(value, 2) for key, value in emp.items()}  # round to 3 decimal places
stats = ''
for k,v in emp.items():
    stats += f'{k}:\n{v}\n'

# plot the breakpoints with probability distribution for each axis
if args.plot:
    g = sns.jointplot(x=x, y=y, kind="scatter", marginal_kws=dict(fill=True), alpha=0.3)
    g.ax_marg_x.get_xaxis().get_offset_text().set_visible(False)
    g.ax_marg_y.get_yaxis().get_offset_text().set_visible(False)
    g.fig.subplots_adjust(right=0.75)  # leave space on the right side
    g.fig.subplots_adjust(left=0.2, bottom=0.2)  # leave space on the left and bottom side

    g.fig.text(0.77, 0.5, stats, va='center', ha='left', fontsize=args.fontsize-5)

    plt.xlabel(f'{args.left} base pair position', fontsize=args.fontsize)
    plt.ylabel(f'{args.right} base pair position', fontsize=args.fontsize)
    # REMOVE the offset text from the marginal axes (top, right)
    plt.savefig(f'{args.output}_scat_dens.png')
    plt.close()


if args.mc:
    # random data
    random.seed(args.seed)

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
        rand = clark_evans_R(df_null['x'].values, df_null['y'].values, x_start=l_start, x_end=l_end, y_start=r_start, y_end=r_end, log10=args.log)
        R_null.append(rand['R'])

    # compute the monte carlo p-value by two-tailed test
    R_null = np.array(R_null)
    mean_null = np.mean(R_null)
    # test whether a null value is as or more extreme than the observed value
    def is_null_extreme(null, observed, mean_null):
        return abs(null - mean_null) >= abs(observed - mean_null)
    numerator = 0
    for x_i in R_null:
        if is_null_extreme(x_i, emp['R'], mean_null):
            numerator += 1
    p_value = (numerator + 1) / (args.permutations + 1)

    if args.plot:
        # do a histogram of null distribution and vline of empirical R stat
        plt.hist(R_null, bins=20, alpha=0.5, color='blue', label='Null distribution')
        plt.axvline(x=emp['R'], color='red', linestyle='dashed', linewidth=2, label=f'Empirical R: {emp["R"]} (p={p_value:.3f})')
        plt.xlabel('R stat', fontsize=args.fontsize)
        plt.ylabel('Frequency', fontsize=args.fontsize)
        plt.title('Null distribution of R stats', fontsize=args.fontsize)
        plt.legend()
        plt.savefig(f'{args.output}_null_dist.png')
