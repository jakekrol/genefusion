#!/usr/bin/env python3
import matplotlib.pyplot as plt
import os,sys
import argparse
import pandas as pd
import math

# input: 1) intersect unswapped file of left gene, 2) names of left and right gene
# output: plot of breakpoints

GENEFILE='/data/jake/genefusion/data/gene_file.txt.latest'
parser = argparse.ArgumentParser(description='Plot breakpoints')
parser.add_argument('-l', '--left', help='left gene name', type=str, required=True)
parser.add_argument('-r', '--right',help='right gene name', type=str, required=True)
parser.add_argument('-i', '--input', help='input intersect file w/ breakpoints', type=str, required=True)
parser.add_argument('-o', '--output', help='output file name', type=str, required=True)
parser.add_argument('-f', '--font_size', help='font size', type=int, default=20)
parser.add_argument('-c', '--color', type=str, default='blue', help='comma separated list of colors')
parser.add_argument('-s', '--size', help='dot size', type=int, default=150)
parser.add_argument('-n', '--n_samples', action='store_true', help='include number of samples in title')
parser.add_argument('-b', '--both', action='store_true', help='add this flag if the input file ' +
    'has an extra col indicating the sample type.') 
args = parser.parse_args()
args.color = args.color.split(',')
assert len(args.color) <=2, f'Only 2 colors are supported. Found {len(args.color)} colors {args.color}'

# lookup start and end positions of genes
df_genes = pd.read_csv(GENEFILE, sep='\t', header=None)
df_genes.columns = ['chrm', 'start', 'end', 'name', 'strand']
l_start = df_genes[df_genes['name'] == args.left]['start'].values[0]
l_end = df_genes[df_genes['name'] == args.left]['end'].values[0]
r_start = df_genes[df_genes['name'] == args.right]['start'].values[0]
r_end = df_genes[df_genes['name'] == args.right]['end'].values[0]

# read
df = pd.read_csv(args.input, sep='\t', header=None)
if not args.both:
    df.columns = [
        'hit_chrm', 'hit_start', 'hit_end', 'hit_name', 'hit_strand',
        'left_chrm', 'left_start', 'left_end', 'left_strand',
        'right_chrm', 'right_start', 'right_end', 'right_strand',
        'xxx', 'sample'
    ]
else:
    df.columns = [
        'hit_chrm', 'hit_start', 'hit_end', 'hit_name', 'hit_strand',
        'left_chrm', 'left_start', 'left_end', 'left_strand',
        'right_chrm', 'right_start', 'right_end', 'right_strand',
        'xxx', 'sample', 'sample_type'
    ]
    df.sort_values(by=['sample_type'], inplace=True)
# filter
df = df[df['hit_name'] == args.right]
m = df.shape[0]
assert m > 0, f'No hits found for {args.right} in {args.input}'
print(f'Found {m} hits for {args.right} in {args.input}')

# get interval midpoints
df['left_mid'] = ((df['left_start'] + df['left_end']) / 2).apply(math.floor)
df['right_mid'] = ((df['right_start'] + df['right_end']) / 2).apply(math.floor)

plt.figure(figsize=(10, 10))
x = df['left_mid'].values
y = df['right_mid'].values
# if args.color:
#     color = df['sample'].values
#     cset = set(color)
#     # plot each sample in different color
#     for c in cset:
#         x1 = x[color == c]
#         y1 = y[color == c]
#         plt.scatter(x1, y1, alpha=0.5, s=args.size, label=c)
#         plt.legend()
if not args.both:
    plt.scatter(x, y, alpha=0.5,s=args.size, color=args.color) 
else:
    n_tumor = df[df['sample_type'] == 'tumor'].shape[0]
    n_normal = df[df['sample_type'] == 'normal'].shape[0]
    t=False
    n=False
    for i,row in df.iterrows():
        if row['sample_type'] == 'tumor':
            if t == False:
                t=True
                plt.scatter(row['left_mid'], row['right_mid'], alpha=0.5, s=args.size, color=args.color[0], label=f'tumor ({n_tumor}/{m})')
            else:
                plt.scatter(row['left_mid'], row['right_mid'], alpha=0.5, s=args.size, color=args.color[0])
        else:
            if n == False:
                n=True
                plt.scatter(row['left_mid'], row['right_mid'], alpha=0.5, s=args.size, color=args.color[1], label=f'normal ({n_normal}/{m})')
            else:
                plt.scatter(row['left_mid'], row['right_mid'], alpha=0.5, s=args.size, color=args.color[1])
    plt.legend()
    
plt.xlim(l_start, l_end)
plt.ylim(r_start, r_end)
plt.xlabel(f'{args.left} breakpoints', fontsize=args.font_size)
plt.ylabel(f'{args.right} breakpoints', fontsize=args.font_size)
plt.title(f'Breakpoint plot for {args.left} and {args.right} ({m} samples)', fontsize=args.font_size)
plt.savefig(args.output, dpi=300)







