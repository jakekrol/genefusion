#!/usr/bin/env python

import gzip
import glob
import pandas as pd
import matplotlib.pyplot as plt

filename_intersect='ERG.giggle.clean.swap.intersect.bed.gz'
intervals_1kg_file = f'./g2f_out/high_coverage_1000g_dna/{filename_intersect}'
intervals_prostate_tumor_file = f'./g2f_out/prostate_tumor/{filename_intersect}'
intervals_prostate_normal_file = f'./g2f_out/prostate_normal/{filename_intersect}'

columns = ['right_gene', 'right_chromosome', 'right_start', 'right_end', 'right_strand', 'left_chromosome', 'left_start', 'left_end', 'left_strand']

df_1kg = pd.read_csv(intervals_1kg_file, sep='\t', header=None, usecols=[3] + list(range(5,13)), comment='#')
df_1kg.columns = columns
df_1kg = df_1kg[df_1kg['right_gene'] == 'TMPRSS2']
df_1kg['strand_config'] = df_1kg.apply(lambda row: f"{row['left_strand']},{row['right_strand']}", axis=1)
df_1kg['cohort'] = '1kg'

df_prostate_tumor = pd.read_csv(intervals_prostate_tumor_file, sep='\t', header=None, usecols= [3] + list(range(5,13)), comment='#')
df_prostate_tumor.columns = columns
df_prostate_tumor = df_prostate_tumor[df_prostate_tumor['right_gene'] == 'TMPRSS2']
df_prostate_tumor['strand_config'] = df_prostate_tumor.apply(lambda row: f"{row['left_strand']},{row['right_strand']}", axis=1)
df_prostate_tumor['cohort'] = 'prostate_tumor'


df_prostate_normal = pd.read_csv(intervals_prostate_normal_file, sep='\t', header=None, usecols= [3] + list(range(5,13)), comment='#')
df_prostate_normal.columns = columns
df_prostate_normal = df_prostate_normal[df_prostate_normal['right_gene'] == 'TMPRSS2']
df_prostate_normal['strand_config'] = df_prostate_normal.apply(lambda row: f"{row['left_strand']},{row['right_strand']}", axis=1)
df_prostate_normal['cohort'] = 'prostate_normal'

df_concat = pd.concat([df_1kg, df_prostate_tumor, df_prostate_normal], ignore_index=True)
# barplot strand configs by cohort
color_map = {
    '1,1': 'lightblue',
    '1,-1': 'red',
    '-1,1': 'green',
    '-1,-1': 'purple'
}
strand_config_counts = df_concat.groupby(['cohort', 'strand_config']).size().reset_index(name='count')
strand_config_counts['color'] = strand_config_counts['strand_config'].map(color_map)
fig, ax = plt.subplots(1,3,figsize=(10, 5))

ax[0].bar(
    strand_config_counts[strand_config_counts['cohort'] == '1kg']['strand_config'],
    strand_config_counts[strand_config_counts['cohort'] == '1kg']['count'],
    color=strand_config_counts[strand_config_counts['cohort'] == '1kg']['color']
)
ax[0].set_title('1kg')
ax[0].set_xlabel('Strand Configuration')
ax[0].set_ylabel('Count')

ax[1].bar(
    strand_config_counts[strand_config_counts['cohort'] == 'prostate_tumor']['strand_config'],
    strand_config_counts[strand_config_counts['cohort'] == 'prostate_tumor']['count'],
    color=strand_config_counts[strand_config_counts['cohort'] == 'prostate_tumor']['color']
)
ax[1].set_title('Prostate Tumor')
ax[1].set_xlabel('Strand Configuration')
ax[1].set_ylabel('Count')

ax[2].bar(
    strand_config_counts[strand_config_counts['cohort'] == 'prostate_normal']['strand_config'],
    strand_config_counts[strand_config_counts['cohort'] == 'prostate_normal']['count'],
    color=strand_config_counts[strand_config_counts['cohort'] == 'prostate_normal']['color']
)
ax[2].set_title('Prostate Normal')
ax[2].set_xlabel('Strand Configuration')
ax[2].set_ylabel('Count')

plt.tight_layout()
plt.savefig('strand_config_counts.png')