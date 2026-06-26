#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import pandas as pd
from polymerization.datasets import *
from polymerization.stix2fusion import *
from polymerization.io import *

parser = argparse.ArgumentParser(description='')
parser.add_argument('--score', default='score.in_tissue.tsv')
parser.add_argument('--outfile', default='hist_score.png')
parser.add_argument('--bed', default='../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed')
parser.add_argument('--bins', type=int, default=15)
args = parser.parse_args()

df_bed = read_bed(args.bed, gene_col_idx=3)
df_score = pd.read_csv(args.score, sep='\t')
df_tumor_fusion =get_pcawg_recurrent_tumor_fusions()
df_tumor_fusion = left_sort_fusion_set(df_tumor_fusion, df_bed=df_bed)
df_score['gene_left'] = df_score['fusion'].apply(lambda x: x.split('--')[0])
df_score['gene_right'] = df_score['fusion'].apply(lambda x: x.split('--')[1])
df_merge = pd.merge(df_score, df_tumor_fusion, how='left', on=['gene_left', 'gene_right'])
df_merge_tumor_reported = df_merge[df_merge['label']==1]
df_merge_tumor_reported = df_merge_tumor_reported[df_merge_tumor_reported['reported_previously']!='False']
df_merge_tumor_unreported = df_merge[df_merge['label']==1]
df_merge_tumor_unreported = df_merge_tumor_unreported[df_merge_tumor_unreported['reported_previously']=='False']
df_merge_normal = df_merge[df_merge['label']==0]
fig, ax = plt.subplots(figsize=(6, 4))

for df in [df_merge_normal, df_merge_tumor_unreported, df_merge_tumor_reported]:
	if df is df_merge_normal:
		label = 'Normal'
		color = 'blue'
	elif df is df_merge_tumor_unreported:
		label = 'Tumor (unreported)'
		color = 'orange'
	else:
		label = 'Tumor (reported)'
		color = 'red'
	bins = int(min(args.bins, len(df['score'].unique())))
	# For very small groups (few points) draw vertical lines
	if len(df['score']) <= 10:
		for score in df['score']:
			ax.axvline(score, color=color, alpha=0.5, label=label)
	else:
		ax.hist(df['score'], bins=bins, alpha=0.5, label=label, color=color)
	print(df['score'].describe())


ax.set_xlabel('Score')
ax.set_ylabel('Count')
fig.legend(loc='right', fontsize=7)
fig.savefig(args.outfile, dpi=300)
df_merge_tumor_unreported.to_csv('score.in_tissue.tumor_unreported.tsv', sep='\t', index=False)

