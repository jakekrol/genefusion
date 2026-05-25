#!/usr/bin/env python

import gzip
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


SEED=0
np.random.seed(SEED)

def get_sample_breakpoints(df):
    sample_ranges = (
        df.groupby('sample').agg(
            left_bp=('left_start', 'min'),
            right_bp=('right_start', 'max')
        )
    )
    sample_ranges.reset_index(inplace=True)
    return sample_ranges

filename_intersect='ERG.giggle.clean.swap.intersect.bed.gz'
intervals_1kg_file = f'./g2f_out/high_coverage_1000g_dna/{filename_intersect}'
intervals_prostate_tumor_file = f'./g2f_out/prostate_tumor/{filename_intersect}'
intervals_prostate_normal_file = f'./g2f_out/prostate_normal/{filename_intersect}'

col_idxs = [3] + list(range(5,13)) + [14]

columns = ['right_gene', 'right_chromosome', 'right_start', 'right_end', 'right_strand', 'left_chromosome', 'left_start', 'left_end', 'left_strand', 'sample']

df_1kg = pd.read_csv(intervals_1kg_file, sep='\t', header=None, usecols=col_idxs, comment='#')
df_1kg.columns = columns
df_1kg = df_1kg[df_1kg['right_gene'] == 'TMPRSS2']
df_1kg['strand_config'] = df_1kg.apply(lambda row: f"{row['left_strand']},{row['right_strand']}", axis=1)
df_1kg['cohort'] = '1kg'
breakpoints_1kg = get_sample_breakpoints(df_1kg)
df_1kg = pd.merge(df_1kg, breakpoints_1kg, on='sample', how='left')

df_prostate_tumor = pd.read_csv(intervals_prostate_tumor_file, sep='\t', header=None, usecols=col_idxs, comment='#')
df_prostate_tumor.columns = columns
df_prostate_tumor = df_prostate_tumor[df_prostate_tumor['right_gene'] == 'TMPRSS2']
df_prostate_tumor['strand_config'] = df_prostate_tumor.apply(lambda row: f"{row['left_strand']},{row['right_strand']}", axis=1)
df_prostate_tumor['cohort'] = 'prostate_tumor'
breakpoints_prostate_tumor = get_sample_breakpoints(df_prostate_tumor)
df_prostate_tumor = pd.merge(df_prostate_tumor, breakpoints_prostate_tumor, on='sample', how='left')


df_prostate_normal = pd.read_csv(intervals_prostate_normal_file, sep='\t', header=None, usecols=col_idxs, comment='#')
df_prostate_normal.columns = columns
df_prostate_normal = df_prostate_normal[df_prostate_normal['right_gene'] == 'TMPRSS2']
df_prostate_normal['strand_config'] = df_prostate_normal.apply(lambda row: f"{row['left_strand']},{row['right_strand']}", axis=1)
df_prostate_normal['cohort'] = 'prostate_normal'
breakpoints_prostate_normal = get_sample_breakpoints(df_prostate_normal)
df_prostate_normal = pd.merge(df_prostate_normal, breakpoints_prostate_normal, on='sample', how='left')

df_concat = pd.concat([df_1kg, df_prostate_tumor, df_prostate_normal], ignore_index=True)

# scatter plot of left and right breakpoints colored by strand config
color_map = {
    '1,1': 'lightblue',
    '1,-1': 'red',
    '-1,1': 'green',
    '-1,-1': 'purple'
}
marker_map = {
	'1,1': '.',
	'1,-1': 'v',
	'-1,1': '^',
	'-1,-1': '<'
}
df_concat['color'] = df_concat['strand_config'].map(color_map)
# add jitter
df_concat['left_bp'] = df_concat['left_bp'] + np.random.normal(0, 1000, size=len(df_concat))
df_concat['right_bp'] = df_concat['right_bp'] + np.random.normal(0, 1000, size=len(df_concat))
cohorts = ['1kg', 'prostate_tumor', 'prostate_normal']
strand_configs = ['1,1', '1,-1', '-1,1', '-1,-1']
fig, ax = plt.subplots(len(cohorts), len(strand_configs), figsize=(14, 9), sharex=True, sharey=True)
for row_index, cohort in enumerate(cohorts):
    cohort_df = df_concat[df_concat['cohort'] == cohort]
    for col_index, strand_config in enumerate(strand_configs):
        axis = ax[row_index, col_index]
        strand_df = cohort_df[cohort_df['strand_config'] == strand_config]
        if not strand_df.empty:
            axis.scatter(
                strand_df['left_bp'],
                strand_df['right_bp'],
                alpha=0.7,
                c=color_map.get(strand_config, 'gray'),
                s=3,
                marker=marker_map.get(strand_config, '.')
            )
        if row_index == 0:
            axis.set_title(strand_config)
        if col_index == 0:
            axis.set_ylabel(f'{cohort}\nTMPRSS2 Breakpoint')
        if row_index == len(cohorts) - 1:
            axis.set_xlabel('ERG Breakpoint')
        axis.tick_params(labelsize=8)
plt.tight_layout()
fig.legend(handles=[plt.Line2D([0], [0], marker='.', color='w', label='1,1', markerfacecolor='lightblue', markersize=5),
					plt.Line2D([0], [0], marker='v', color='w', label='1,-1', markerfacecolor='red', markersize=5),
					plt.Line2D([0], [0], marker='^', color='w', label='-1,1', markerfacecolor='green', markersize=5),
					plt.Line2D([0], [0], marker='<', color='w', label='-1,-1', markerfacecolor='purple', markersize=5)], loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
plt.savefig('breakpoint_scatter_strand_config.png')
