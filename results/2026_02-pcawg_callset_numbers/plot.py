#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
import os,sys

parser=argparse.ArgumentParser(description='Plot PCAWG fusion callset stats')
parser.add_argument('-i', '--input', help='pcawg fusion callset tsv file', 
                    default='../../data/2025_04-pcawg_fusions/gene.fusions.V1.tsv')
parser.add_argument('-o', '--outdir', help='output directory for plots', required=True)
args=parser.parse_args()

df = pd.read_csv(args.input, sep='\t')

### plot num. unique fusions per histology
df['histology_simple'] = df['histology_abbreviation'].apply(lambda x: x.split('-')[0].lower())
histology_cts = df.groupby('histology_simple')['fusion_id'].nunique().reset_index()
# bar plot it sort by num. fusions and color it black
# annotate the total number of fusions in the upper right region
# display number at top of each bar
total=histology_cts['fusion_id'].sum()
histology_cts = histology_cts.sort_values('fusion_id', ascending=False)
plt.figure(figsize=(10,6))
bars = plt.bar(histology_cts['histology_simple'], histology_cts['fusion_id'], color='black')
plt.annotate(f'Total: {total}', xy=(0.9, 0.9), xycoords='axes fraction', ha='right', va='top', fontsize=12)
# add number labels on top of each bar
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height)}',
             ha='center', va='bottom')
plt.xticks(rotation=90)
plt.xlabel('Histology')
plt.ylabel('Number of unique fusions')
plt.title('PCAWG fusion callset')
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, 'num_fusions_per_histology.png'))

### Bar plot of sample counts, log scale y-axis
# horizontal axis is num samples with fusion, vertical axis is num fusions with that many samples
# ticks should be exhaustive
# count number of samples per fusion
sample_counts = df.groupby('fusion_id')['icgc_donor_id'].nunique().reset_index()
# count number of fusions per sample count
sample_counts = sample_counts['icgc_donor_id'].value_counts().reset_index()
sample_counts.columns = ['num_samples', 'num_fusions']
plt.figure(figsize=(10,6))
ax = plt.gca()
plt.bar(sample_counts['num_samples'], sample_counts['num_fusions'], color='black')
ax.set_xticks(sample_counts['num_samples'])
plt.yscale('log')
plt.xlabel('Number of samples with fusion')
plt.ylabel('Number of fusions')
plt.title('PCAWG fusion callset sample counts')
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, 'samples_per_fusion.png'))