#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
FONTSIZE=9
SCORE_COL='score-w_1000g-wnormal_0.5'
f_evidence='score-all-tumor-normal-recurrent.tsv'
# f_roc = './roc_data-neg_normal/roc_input-neg_normal-w_0.5-w_1000g.tsv'
# df_roc=pd.read_csv(f_roc,sep='\t')
df=pd.read_csv(f_evidence,sep='\t')
# df = pd.merge(df_evidence, df_roc, on=['gene_left', 'gene_right'])

# see f_roc for how these were selected (low scores)
# tail -n +2 ./roc_data-neg_normal/roc_input-neg_normal-w_0.5-w_1000g.tsv | sort -k3 -gr
low_score_tumor_fusions=set(
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
df['is_low_score_tumor_fusion'] = 0
for i,row in df.iterrows():
    gl = row['gene_left']
    gr = row['gene_right']
    x = (gl, gr)
    if x in low_score_tumor_fusions:
        indices.append(i)
        df.loc[i,'is_low_score_tumor_fusion'] = 1
# summarize tumor,normal read,sample evidence for each fusion
# onekg_cols=[46,47,28,29]
onekg_cols = df.columns[df.columns.str.contains('1000g').tolist()].tolist()
# drop score columns from onekg_cols
onekg_cols = [x for x in onekg_cols if 'score' not in x]
onekg_reads_cols = [x for x in onekg_cols if 'reads' in x]
onekg_samples_cols = [x for x in onekg_cols if 'samples' in x]
normal_cols=df.columns[df.columns.str.contains('normal').tolist()].tolist()
tumor_cols=df.columns[df.columns.str.contains('tumor').tolist()].tolist()
normal_cols_reads = [x for x in normal_cols + onekg_cols if 'reads' in x]
normal_cols_reads_no_1000g = [x for x in normal_cols_reads if '1000g' not in x]
normal_cols_samples = [x for x in normal_cols + onekg_cols if 'samples' in x]
normal_cols_samples_no_1000g = [x for x in normal_cols_samples if '1000g' not in x]
tumor_cols_reads = [x for x in tumor_cols if 'reads' in x]
tumor_cols_samples = [x for x in tumor_cols if 'samples' in x]
df['total_normal_reads'] = df[normal_cols_reads].sum(axis=1).values
df['total_normal_samples'] = df[normal_cols_samples].sum(axis=1).values
df['total_normal_reads_no_1000g'] = df[normal_cols_reads_no_1000g].sum(axis=1).values
df['total_normal_samples_no_1000g'] = df[normal_cols_samples_no_1000g].sum(axis=1).values
df['total_tumor_reads'] = df[tumor_cols_reads].sum(axis=1).values
df['total_tumor_samples'] = df[tumor_cols_samples].sum(axis=1).values
df['total_1000g_reads'] = df[onekg_reads_cols].sum(axis=1).values
df['total_1000g_samples'] = df[onekg_samples_cols].sum(axis=1).values

evidence_cols = [
    'total_normal_reads',
    'total_normal_reads_no_1000g',
    'total_1000g_reads',
    'total_tumor_reads',
    'total_normal_samples',
    'total_normal_samples_no_1000g',
    'total_1000g_samples',
    'total_tumor_samples'
]

fig, ax = plt.subplots(2,4, figsize=(10,5))
for i, col in enumerate(evidence_cols):
    row = i // 4
    col_idx = i % 4
    def color_scheme(label, is_low_score):
        if label == 0:
            return 'tab:blue'
        elif label == 1:
            if is_low_score == 1:
                return 'tab:orange'
            elif is_low_score == 0:
                return 'tab:red'
            else:
                raise ValueError(f'Unexpected is_low_score value: {is_low_score}')
        else:
            raise ValueError(f'Unexpected label: {label}')
    df['color'] = df.apply(lambda row: color_scheme(row['label'], row['is_low_score_tumor_fusion']), axis=1).values
    
    ax[row, col_idx].scatter(df[col], df[SCORE_COL], c=df['color'], s=3)
    
    ax[row, col_idx].set_ylabel("Score", fontsize=FONTSIZE)
    # if reads log scale
    if 'reads' in col:
        ax[row, col_idx].set_xscale('log')
        ax[row, col_idx].set_xlabel(f'Log scale {col.replace("_", " ")}', fontsize=FONTSIZE)
    else:
        ax[row, col_idx].set_xlabel(col.replace('_', ' ').capitalize(), fontsize=FONTSIZE)
    low_read_depth_max = df[df['is_low_score_tumor_fusion'] == 1][col].max()
    # rm top and right spines
    ax[row, col_idx].spines['top'].set_visible(False)
    ax[row, col_idx].spines['right'].set_visible(False)
# legend describing colors
legend_handles = [
    Line2D([0], [0], marker='o', color='w', label='Normal fusion', markerfacecolor='tab:blue', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='High score tumor fusion', markerfacecolor='tab:red', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Low score tumor fusion', markerfacecolor='tab:orange', markersize=8),
]
plt.tight_layout(rect=[0, 0, 1, 0.85])
fig.legend(handles=legend_handles, title='', frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.98))
plt.savefig('scatter-evidence_score.png', dpi=300)

# save relevant plot data
plot_data_cols = ['gene_left', 'gene_right', SCORE_COL] + evidence_cols + ['label', 'is_low_score_tumor_fusion']
df[plot_data_cols].to_csv('scatter-evidence_score_data.tsv', sep='\t', index=False)
# save detailed plot data
raw_plot_data_cols = ['gene_left', 'gene_right', SCORE_COL] + df.columns[df.columns.str.contains('reads|samples')].tolist() + ['label', 'is_low_score_tumor_fusion']
df[raw_plot_data_cols].to_csv('scatter-evidence_score_raw_data.tsv', sep='\t', index=False)



# df_sub = df[df['label_y'] == 1]
# score_col='score-w_1000g-wnormal_0.5_y'
# fig, ax = plt.subplots()
# x=df_sub['total_tumor_reads']
# y=df_sub[score_col]
# color=df_sub['is_low_score_tumor_fusion'].map({0: 'tab:blue', 1: 'tab:orange'}).values
# ax.scatter(x,y, c=color)
# ax.set_xlabel('Total tumor read depth', fontsize=FONTSIZE)
# ax.set_ylabel("Score", fontsize=FONTSIZE)
# ax.set_xscale('log')
# low_read_depth_max = df_sub[df_sub['is_low_score_tumor_fusion'] == 1]['total_tumor_reads'].max()
# ax.axvline(low_read_depth_max, color='gray', linestyle='--', label=f'Max tumor read depth: {low_read_depth_max:.0f}')
# # rm top and right spines

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# # legend describing colors
# legend_handles = [
#     Line2D([0], [0], marker='o', color='w', label='High scoring tumor fusion', markerfacecolor='tab:blue', markersize=8),
#     Line2D([0], [0], marker='o', color='w', label='Low scoring tumor fusion', markerfacecolor='tab:orange', markersize=8),
# ]
# ax.legend(handles=legend_handles, title='', frameon=False, loc='upper center')


# plt.savefig('scatter-tumor_reads_score.png')

# fig, ax = plt.subplots()
# #x=np.log10(df['total_normal_reads'])
# #y=np.log10(df['total_tumor_reads'])
# x = df['total_normal_reads']
# y = df['total_tumor_reads']
# m = max(x.max(), y.max())
# color=df['is_low_score_tumor_fusion'].values
# ax.scatter(x,y, c=color)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('total normal reads')
# ax.set_ylabel('total tumor reads')
# ax.plot(np.linspace(0,m), np.linspace(0,m))

# plt.savefig('scatter-total_reads.png')
