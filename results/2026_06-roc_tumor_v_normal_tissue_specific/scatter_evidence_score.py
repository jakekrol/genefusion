#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from polymerization.stix2fusion import left_sort_fusion_set
from polymerization.io import read_bed

parser = argparse.ArgumentParser(description='Scatter plot of evidence v score')
parser.add_argument('--data', default='all_evidence_with_score.tsv')
parser.add_argument('--tumor_tissues', default='../2026_04-s2f_pcawg_recurrent/recurrent_tumor_fusions.tsv')
parser.add_argument('--normal_tissues', default='../2026_06-g2f-normal_fusions_tissue_specific/normal_fusion_tissue_map.tsv')
parser.add_argument('--tumor_bed', default='../2026_05-g2f-tumor_and_normal-rerun/grch37.genes.bed.added')
parser.add_argument('--output_scatter', default='evidence_score_scatter.png')
parser.add_argument('--output_data', default='evidence_score_scatter_data.tsv')
parser.add_argument('--marker_size', type=float, default=5, help='marker size for scatter plot')
args = parser.parse_args()

df = pd.read_csv(args.data, sep='\t')
df.drop(columns=['reads_low_coverage_1000g_dna', 'samples_low_coverage_1000g_dna'], inplace=True)
df_tumor_bed = read_bed(args.tumor_bed, gene_col_idx=3)
df_tumor_tissues = pd.read_csv(args.tumor_tissues, sep='\t')
df_tumor_tissues = left_sort_fusion_set(df_tumor_tissues, df_tumor_bed)
df_normal_tissues = pd.read_csv(args.normal_tissues, sep='\t')
def rename_tissues(str_tissues):
    tissues = []
    for t in str_tissues.split(","):
        if t == 'gall bladder':
            tissues.append('gallbladder')
        elif t == 'bone_marrow':
            tissues.append('bone')
        else:
            tissues.append(t)
    return ",".join(tissues)
df_normal_tissues['tissues'] = df_normal_tissues['tissues'].apply(rename_tissues)
df_tissues = pd.concat([df_tumor_tissues[['gene_left', 'gene_right', 'tissues']], df_normal_tissues], ignore_index=True)
df = pd.merge(
    df,
    df_tissues,
    on=['gene_left', 'gene_right'],
    how='left'
)
df.drop_duplicates(subset=['gene_left', 'gene_right', 'label'], inplace=True)
df = df.sort_values(by='score', ascending=False)

cols = df.columns.tolist()
# read evidence columns
read_cols = [c for c in cols if c.startswith('reads_')]
tumor_read_cols = [c for c in read_cols if 'tumor' in c]
normal_read_cols = [c for c in read_cols if 'normal' in c]
normal_read_cols = normal_read_cols +  [c for c in read_cols if '1000g' in c]
# sample evidence columns
sample_cols = [c for c in cols if c.startswith('samples_')]
tumor_sample_cols = [c for c in sample_cols if 'tumor' in c]
normal_sample_cols = [c for c in sample_cols if 'normal' in c]
normal_sample_cols = normal_sample_cols + [c for c in sample_cols if '1000g' in c]
def fusion2evidence(row, evidence_cols, normal=False):
    # get tissues
    tissues = row['tissues'].split(",") if pd.notna(row['tissues']) else []
    # subset columns
    row = row[evidence_cols]
    evidence = 0
    # sum evidence across relevant columns
    for idx, val in row.items():
        # for normal evidence always include 1000g
        if '1000g' in idx and normal:
            evidence += val
            continue
        idx_tissue = idx.split("_")[1]
        for t in tissues:
            if t == idx_tissue:
                evidence += val
    return evidence
def fusion2evidence_data(row, evidence_cols):
    tissues = row['tissues'].split(",") if pd.notna(row['tissues']) else []
    data = dict()
    for idx, val in row.items():
        # always include 1000g
        if '1000g' in idx:
            data[idx] = val
            continue
        if idx in evidence_cols:
            idx_tissue = idx.split("_")[1]
            for t in tissues:
                if t == idx_tissue:
                    data[idx] = val
    return data
df['tumor_read_evidence'] = df.apply(lambda row: fusion2evidence(row, tumor_read_cols, normal=False), axis=1)
df['normal_read_evidence'] = df.apply(lambda row: fusion2evidence(row, normal_read_cols, normal=True), axis=1)
df['tumor_sample_evidence'] = df.apply(lambda row: fusion2evidence(row, tumor_sample_cols, normal=False), axis=1)
df['normal_sample_evidence'] = df.apply(lambda row: fusion2evidence(row, normal_sample_cols, normal=True), axis=1)

fig, ax = plt.subplots(2,2, figsize=(10,10))

ax[0,0].scatter(np.log10(df['tumor_read_evidence'] + 1), df['score'], c=df['label'].map({0:'blue', 1:'red'}), alpha=0.5, s=args.marker_size)
ax[0,0].set_xlabel('log10(Tumor Read Evidence)')
ax[0,0].set_ylabel('Score')
ax[0,0].set_title('Tumor Read Evidence vs Score')
ax[0,1].scatter(np.log10(df['normal_read_evidence'] + 1), df['score'], c=df['label'].map({0:'blue', 1:'red'}), alpha=0.5, s=args.marker_size)
ax[0,1].set_xlabel('log10(Normal Read Evidence)')
ax[0,1].set_ylabel('Score')
ax[0,1].set_title('Normal Read Evidence vs Score')
ax[1,0].scatter(df['tumor_sample_evidence'], df['score'], c=df['label'].map({0:'blue', 1:'red'}), alpha=0.5, s=args.marker_size)
ax[1,0].set_xlabel('Tumor Sample Evidence')
ax[1,0].set_ylabel('Score')
ax[1,0].set_title('Tumor Sample Evidence vs Score')
ax[1,1].scatter(df['normal_sample_evidence'], df['score'], c=df['label'].map({0:'blue', 1:'red'}), alpha=0.5, s=args.marker_size)
ax[1,1].set_xlabel('Normal Sample Evidence')
ax[1,1].set_ylabel('Score')
ax[1,1].set_title('Normal Sample Evidence vs Score')

plt.tight_layout()
plt.savefig(args.output_scatter)

df['data'] = df.apply(lambda row: fusion2evidence_data(row, tumor_read_cols + normal_read_cols + tumor_sample_cols + normal_sample_cols), axis=1)   
df.to_csv(args.output_data, sep='\t', index=False)
