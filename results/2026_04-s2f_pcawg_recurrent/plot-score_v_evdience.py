#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
FONTSIZE=15
f_evidence='../2026_04-roc_tumor_and_normal_recurrent_fusions/score-all-tumor-normal-recurrent.tsv'
f_roc = '../2026_04-roc_tumor_and_normal_recurrent_fusions/roc_data-neg_normal/roc_input-neg_normal-w_0.5-w_1000g.tsv'
df_roc=pd.read_csv(f_roc,sep='\t')
df_evidence=pd.read_csv(f_evidence,sep='\t')
df = pd.merge(df_evidence, df_roc, on=['gene_left', 'gene_right'])
fusions=set(
    [
        ('GNS', 'NUP107'),
        ('ACTG2','DGUOK'),
        ('PMEPA1', 'MSL2'),
        ('FBF1', 'AFM'),
        ('TRAF3IP2', 'FYN'),
        ('TACC3', 'FGFR3')
    ]
)
indices = []
df['indicator'] = 0
for i,row in df.iterrows():
    gl = row['gene_left']
    gr = row['gene_right']
    x = (gl, gr)
    if x in fusions:
        indices.append(i)
        df.loc[i,'indicator'] = 1
onekg_cols=[46,47,28,29]
onekg_cols=df.columns[onekg_cols].tolist()
normal_cols=df.columns[df.columns.str.contains('normal').tolist()].tolist()
tumor_cols=df.columns[df.columns.str.contains('tumor').tolist()].tolist()
normal_cols_reads = [x for x in normal_cols + onekg_cols if 'reads' in x]
normal_cols_samples = [x for x in normal_cols + onekg_cols if 'samples' in x]
tumor_cols_reads = [x for x in tumor_cols if 'reads' in x]
tumor_cols_samples = [x for x in tumor_cols if 'sample' in x]
df['total_normal_reads'] = df[normal_cols_reads].sum(axis=1).values
df['total_normal_samples'] = df[normal_cols_samples].sum(axis=1).values
df['total_tumor_reads'] = df[tumor_cols_reads].sum(axis=1).values
df['total_tumor_samples'] = df[tumor_cols_samples].sum(axis=1).values


df_sub = df[df['label_y'] == 1]
score_col='score-w_1000g-wnormal_0.5_y'
fig, ax = plt.subplots()
x=df_sub['total_tumor_reads']
y=df_sub[score_col]
color=df_sub['indicator'].map({0: 'tab:blue', 1: 'tab:orange'}).values
ax.scatter(x,y, c=color)
ax.set_xlabel('Total tumor read depth', fontsize=FONTSIZE)
ax.set_ylabel("Score", fontsize=FONTSIZE)
ax.set_xscale('log')
low_read_depth_max = df_sub[df_sub['indicator'] == 1]['total_tumor_reads'].max()
ax.axvline(low_read_depth_max, color='gray', linestyle='--', label=f'Max tumor read depth: {low_read_depth_max:.0f}')
# rm top and right spines

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# legend describing colors
legend_handles = [
    Line2D([0], [0], marker='o', color='w', label='High scoring tumor fusion', markerfacecolor='tab:blue', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Low scoring tumor fusion', markerfacecolor='tab:orange', markersize=8),
]
ax.legend(handles=legend_handles, title='', frameon=False, loc='upper center')


plt.savefig('scatter-tumor_reads_score.png')

# fig, ax = plt.subplots()
# #x=np.log10(df['total_normal_reads'])
# #y=np.log10(df['total_tumor_reads'])
# x = df['total_normal_reads']
# y = df['total_tumor_reads']
# m = max(x.max(), y.max())
# color=df['indicator'].values
# ax.scatter(x,y, c=color)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('total normal reads')
# ax.set_ylabel('total tumor reads')
# ax.plot(np.linspace(0,m), np.linspace(0,m))

# plt.savefig('scatter-total_reads.png')
