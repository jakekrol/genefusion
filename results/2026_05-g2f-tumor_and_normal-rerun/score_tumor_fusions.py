#!/usr/bin/env python3

import argparse
import pandas as pd
from polymerization.score import *
from polymerization.io import *
import yaml

parser = argparse.ArgumentParser(description='Calculate fusion scores')
parser.add_argument('--fusion_evidence', default='query-evidence-filled.tsv')
parser.add_argument('--tumor_fusions', default='recurrent_tumor_fusions.sort.labeled.tsv')
parser.add_argument('--colmap', default='score_column_map.yaml')
parser.add_argument('--output', default='tumor-scored.tsv')
parser.add_argument('--tumor_fusion_tissues', default='../2026_04-s2f_pcawg_recurrent/recurrent_tumor_fusions.tsv')
parser.add_argument('--shardfile', default='./shardfile.tsv')

args = parser.parse_args()

def tissue_w_data(tissues):
	t_w_data = []
	for t in tissues.split(","):
		for c in categories:
			category_tissue = c.split("_")[0].lower()
			# string check
			if t == category_tissue:
				t_w_data.append(t)
				break
	return t_w_data

# readinputs
df_evidence = pd.read_csv(args.fusion_evidence, sep='\t')
drop_cols = ['reads_low_coverage_1000g_dna', 'samples_low_coverage_1000g_dna']
df_evidence = df_evidence.drop(columns=drop_cols)
df_tumor = pd.read_csv(args.tumor_fusions, sep='\t')
df_tumor_fusion_tissues = pd.read_csv(args.tumor_fusion_tissues, sep='\t')
df_shard = read_giggle_shardfile(args.shardfile)

# get fusions with tissue data
categories = set(df_shard["category"])
df_tumor_fusion_tissues['fusion_key'] = df_tumor_fusion_tissues.apply(
    lambda row: "".join(sorted((row['gene_x'] + row['gene_y']).upper())),
    axis=1,
)
df_tumor['fusion_key'] = df_tumor.apply(
    lambda row: "".join(sorted((row['gene_left'] + row['gene_right']).upper())),
    axis=1,
)
df_tumor = pd.merge(
    df_tumor,
    df_tumor_fusion_tissues[['fusion_key', 'tissues']],
    on='fusion_key',
    how='left'
)
df_tumor['tissue_w_data'] = df_tumor['tissues'].apply(tissue_w_data)
df_tumor = df_tumor[df_tumor['tissue_w_data'].apply(len) > 0].reset_index(drop=True)


# normalize evidence
with open(args.colmap) as f:
    column_map = yaml.safe_load(f)
df_evidence_tumor = pd.merge(
    df_tumor,
    df_evidence,
    on=['gene_left', 'gene_right'],
    how='left'
)
df_evidence_tumor_norm = normalize_evidence_columns(df_evidence_tumor, column_map)
# scale
evidence_cols = []
for col in df_evidence_tumor_norm.columns:
    if (col.startswith('reads')) or (col.startswith('samples')):
        evidence_cols.append(col)
for col in evidence_cols:
    df_evidence_tumor_norm[col] = df_evidence_tumor_norm[col] * 0.5

# score in 3 ways
# i) agg all weighted by subpop size
# ii) agg all uniformly wighted subpop sizes
# iii) score in-tissue only

evidence_categories = set()
for key in column_map.keys():
    ec = "_".join(key.split("_")[1:])
    evidence_categories.add(ec)
### i)
total_tumor_samples = 0
total_normal_samples = 0
tumor_sample_cols = []
normal_sample_cols = []
for key in column_map.keys():
    # only loop over samples column to avoid double counting
    if 'samples' in key:
        if column_map[key]['specimen'] == 'tumor':
            tumor_sample_cols.append(key)
            total_tumor_samples += column_map[key]['total_samples']
        if column_map[key]['specimen'] == 'normal':
            normal_sample_cols.append(key)
            total_normal_samples += column_map[key]['total_samples']
# weigh by subpop size
for key_column_in, sub_dict in column_map.items():
    specimen = sub_dict['specimen']
    total_samples = sub_dict['total_samples']
    if specimen == 'tumor':
        df_evidence_tumor_norm[key_column_in] = df_evidence_tumor_norm[key_column_in] * (total_samples / (total_tumor_samples))
    if specimen == 'normal':
        df_evidence_tumor_norm[key_column_in] = df_evidence_tumor_norm[key_column_in] * (total_samples / (total_normal_samples))

# normal columns get negative weight
for col in df_evidence_tumor_norm.columns:
    if 'normal' in col:
        df_evidence_tumor_norm[col] = df_evidence_tumor_norm[col] * -1
df_evidence_tumor_norm['score_weighted'] = df_evidence_tumor_norm[column_map.keys()].sum(axis=1)
df_evidence_tumor_norm = df_evidence_tumor_norm.sort_values(by='score_weighted', ascending=False).reset_index(drop=True)
df_evidence_tumor_norm['fusion'] = df_evidence_tumor_norm['gene_left'] + '--' + df_evidence_tumor_norm['gene_right']
df_evidence_tumor_norm[['fusion', 'score_weighted']].to_csv('tumor_score.agg.subpop.weighted.tsv', sep='\t', index=False)

### ii)
# reset values
df_evidence_tumor_norm = normalize_evidence_columns(df_evidence_tumor, column_map)
# scale
for col in evidence_cols:
    df_evidence_tumor_norm[col] = df_evidence_tumor_norm[col] * 0.5
# weight
n_tumor_subpops = len(tumor_sample_cols)
n_normal_subpops = len(normal_sample_cols)
for key_column_in, sub_dict in column_map.items():
    specimen = sub_dict['specimen']
    if specimen == 'tumor':
        df_evidence_tumor_norm[key_column_in] = df_evidence_tumor_norm[key_column_in] * (1 / n_tumor_subpops)
    if specimen == 'normal':
        df_evidence_tumor_norm[key_column_in] = df_evidence_tumor_norm[key_column_in] * (1 / n_normal_subpops)
# normal columns get negative weight
for col in df_evidence_tumor_norm.columns:
    if 'normal' in col:
        df_evidence_tumor_norm[col] = df_evidence_tumor_norm[col] * -1
df_evidence_tumor_norm['score_uniform'] = df_evidence_tumor_norm[column_map.keys()].sum(axis=1)
df_evidence_tumor_norm = df_evidence_tumor_norm.sort_values(by='score_uniform', ascending=False).reset_index(drop=True)
df_evidence_tumor_norm['fusion'] = df_evidence_tumor_norm['gene_left'] + '--' + df_evidence_tumor_norm['gene_right']
df_evidence_tumor_norm[['fusion', 'score_uniform']].to_csv('tumor_score.agg.subpop.uniform.tsv', sep='\t', index=False)

### iii)
# reset
df_evidence_tumor_norm = normalize_evidence_columns(df_evidence_tumor, column_map)
# scale
for col in evidence_cols:
    df_evidence_tumor_norm[col] = df_evidence_tumor_norm[col] * 0.5
# for each fusion select only the relevance evidence types from the column maps
for i, row in df_evidence_tumor_norm.iterrows():
    tissues = row['tissue_w_data']
    evidence_cols_tumor = []
    evidence_cols_normal = []
    n_tumor_subpops = 0
    n_normal_subpops = 0
    for col in column_map.keys():
        # manually add 1000g normal evidence if relevant
        if '1000g' in col:
            evidence_cols_normal.append(col)
            n_normal_subpops += 1
            continue
        tissue_key = col.split("_")[1].lower()
        if tissue_key in tissues:
            if 'tumor' in col:
                evidence_cols_tumor.append(col)
                n_tumor_subpops += 1
            if 'normal' in col:
                evidence_cols_normal.append(col)
                n_normal_subpops += 1
    # weight uniformly by number of relevant subpops
    for col in evidence_cols_tumor:
        df_evidence_tumor_norm.at[i, col] = df_evidence_tumor_norm.at[i, col] * (1 / n_tumor_subpops)
    # normal columns get negative weight
    for col in evidence_cols_normal:
        df_evidence_tumor_norm.at[i, col] = df_evidence_tumor_norm.at[i, col] * (1 / n_normal_subpops) * -1
    df_evidence_tumor_norm.at[i, 'score_in_tissue'] = df_evidence_tumor_norm.loc[i, evidence_cols_tumor + evidence_cols_normal].sum()
df_evidence_tumor_norm = df_evidence_tumor_norm.sort_values(by='score_in_tissue', ascending=False).reset_index(drop=True)
df_evidence_tumor_norm['fusion'] = df_evidence_tumor_norm['gene_left'] + '--' + df_evidence_tumor_norm['gene_right']
df_evidence_tumor_norm[['fusion', 'score_in_tissue']].to_csv('tumor_score.in_tissue.tsv', sep='\t', index=False)


        



    


