#!/usr/bin/env python

import numpy as np
import pandas as pd
from polymerization.io import *
from polymerization.analysis import *
import seaborn as sns
import matplotlib.pyplot as plt

BINS_PER_DIM = 10
df_bed = read_bed('../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed', gene_col_idx=3)
ERGSTART=df_bed[df_bed['gene_name'] == 'ERG']['start'].values[0]
ERGEND=df_bed[df_bed['gene_name'] == 'ERG']['end'].values[0]
GENERIGHT='TMPRSS2'
TMPRSS2START=df_bed[df_bed['gene_name'] == GENERIGHT]['start'].values[0]
TMPRSS2END=df_bed[df_bed['gene_name'] == GENERIGHT]['end'].values[0]
_, df_tumor = read_g2f_intersect('../2026_05-erg_tmprss2-strands/g2f_out/prostate_tumor/ERG.giggle.clean.swap.intersect.bed.gz', bgzip=True)
_, df_normal = read_g2f_intersect('../2026_05-erg_tmprss2-strands/g2f_out/prostate_normal/ERG.giggle.clean.swap.intersect.bed.gz', bgzip=True)
_, df_thousg = read_g2f_intersect('../2026_05-erg_tmprss2-strands/g2f_out/high_coverage_1000g_dna/ERG.giggle.clean.swap.intersect.bed.gz', bgzip=True)

# may need to group by sample later to get sample count
# trying read count first
df_tumor_bp = intersect2breakpoints(df_tumor, GENERIGHT)
df_normal_bp = intersect2breakpoints(df_normal, GENERIGHT)
df_thousg_bp = intersect2breakpoints(df_thousg, GENERIGHT)

def bin_breakpoints(df_bp, x_edges, y_edges):
	'''
	bin breakpoints into 2D histogram
	'''
	hist, _, _ = np.histogram2d(df_bp['left_breakpoint'], df_bp['right_breakpoint'], bins=(x_edges, y_edges))
	# swap axes to have ERG on x and TMPRSS2 on y
	# np.histogram2d array puts x on rows and y on columns, so need to transpose to have ERG on x and TMPRSS2 on y
	hist = hist.T
	return hist

# get tumor bins
x_edges = np.linspace(ERGSTART, ERGEND, BINS_PER_DIM + 1)
y_edges = np.linspace(TMPRSS2START, TMPRSS2END, BINS_PER_DIM + 1)
tumor_hist = bin_breakpoints(df_tumor_bp, x_edges, y_edges)
# get normal bins
normal_hist = bin_breakpoints(df_normal_bp, x_edges, y_edges)
# get thousg bins
thousg_hist = bin_breakpoints(df_thousg_bp, x_edges, y_edges)
# plot all 3 histograms in 1 row, 3 columns
xticklabels = np.round(x_edges[:-1] * 1e-6, 2)
yticklabels = np.round(y_edges[:-1] * 1e-6, 2)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Tumor
sns.heatmap(tumor_hist, xticklabels=xticklabels, yticklabels=yticklabels, cmap='Reds', 
            cbar_kws={'label': 'Read Count'}, ax=axes[0])
axes[0].set_xlabel('ERG Breakpoint (Mb)')
axes[0].set_ylabel('TMPRSS2 Breakpoint (Mb)')
axes[0].set_title('Tumor')

# Normal
sns.heatmap(normal_hist, xticklabels=xticklabels, yticklabels=yticklabels, cmap='Blues', 
            cbar_kws={'label': 'Read Count'}, ax=axes[1])
axes[1].set_xlabel('ERG Breakpoint (Mb)')
axes[1].set_ylabel('TMPRSS2 Breakpoint (Mb)')
axes[1].set_title('Normal')

# 1000G
sns.heatmap(thousg_hist, xticklabels=xticklabels, yticklabels=yticklabels, cmap='Greens', 
            cbar_kws={'label': 'Read Count'}, ax=axes[2])
axes[2].set_xlabel('ERG Breakpoint (Mb)')
axes[2].set_ylabel('TMPRSS2 Breakpoint (Mb)')
axes[2].set_title('1000G')

plt.tight_layout()
plt.savefig('tumor_breakpoint_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

# breakpoint aware normal reads
N = normal_hist + thousg_hist
N_bp_aware = np.minimum(N, tumor_hist)
num_N_bp = N_bp_aware.sum()

# number of non-breakpoint aware normal reads 
num_N = N.sum()

# barplot of breakpoint aware vs non-breakpoint aware normal reads
labels = ['Normal reads', 'Normal reads breakpoint aware']
counts = [num_N, num_N_bp]
plt.figure(figsize=(6, 4))
sns.barplot(x=labels, y=counts, palette=['blue', 'orange'])
plt.ylabel('Read Count')
plt.title('Breakpoint Aware vs Non-Breakpoint Aware Normal Reads')
plt.savefig('normal_breakpoint_aware_barplot.png', dpi=300, bbox_inches='tight')
plt.close()

# now test sample count

df_tumor_bp_sample = intersect2breakpoints(df_tumor, GENERIGHT, group_by_sample=True)
df_normal_bp_sample = intersect2breakpoints(df_normal, GENERIGHT, group_by_sample=True)
df_1000g_bp_sample = intersect2breakpoints(df_thousg, GENERIGHT, group_by_sample=True)

# repeat analysis
tumor_hist_sample = bin_breakpoints(df_tumor_bp_sample, x_edges, y_edges)
normal_hist_sample = bin_breakpoints(df_normal_bp_sample, x_edges, y_edges)
thousg_hist_sample = bin_breakpoints(df_1000g_bp_sample, x_edges, y_edges)

# plot all 3 histograms in 1 row, 3 columns
xticklabels = np.round(x_edges[:-1] * 1e-6, 2)
yticklabels = np.round(y_edges[:-1] * 1e-6, 2)
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
# Tumor
sns.heatmap(tumor_hist_sample, xticklabels=xticklabels, yticklabels=yticklabels, cmap='Reds', 
			cbar_kws={'label': 'Sample Count'}, ax=axes[0])
axes[0].set_xlabel('ERG Breakpoint (Mb)')
axes[0].set_ylabel('TMPRSS2 Breakpoint (Mb)')
axes[0].set_title('Tumor')
# Normal
sns.heatmap(normal_hist_sample, xticklabels=xticklabels, yticklabels=yticklabels, cmap='Blues', 
			cbar_kws={'label': 'Sample Count'}, ax=axes[1])
axes[1].set_xlabel('ERG Breakpoint (Mb)')
axes[1].set_ylabel('TMPRSS2 Breakpoint (Mb)')
axes[1].set_title('Normal')
# 1000G
sns.heatmap(thousg_hist_sample, xticklabels=xticklabels, yticklabels=yticklabels, cmap='Greens', 
			cbar_kws={'label': 'Sample Count'}, ax=axes[2])
axes[2].set_xlabel('ERG Breakpoint (Mb)')
axes[2].set_ylabel('TMPRSS2 Breakpoint (Mb)')
axes[2].set_title('1000G')
plt.tight_layout()
plt.savefig('tumor_breakpoint_heatmap_sample.png', dpi=300, bbox_inches='tight')
plt.close()

# compare counts
N_sample = normal_hist_sample + thousg_hist_sample
N_sample_bp_aware = np.minimum(N_sample, tumor_hist_sample)
num_N_sample_bp = N_sample_bp_aware.sum()
num_N_sample = N_sample.sum()
labels = ['Normal samples', 'Normal samples breakpoint aware']
counts = [num_N_sample, num_N_sample_bp]
plt.figure(figsize=(6, 4))
sns.barplot(x=labels, y=counts, palette=['blue', 'orange'])
plt.ylabel('Sample Count')
plt.title('Breakpoint Aware vs Non-Breakpoint Aware Normal Samples')
plt.savefig('normal_sample_breakpoint_aware_barplot.png', dpi=300, bbox_inches='tight')
plt.close()



breakpoint()





