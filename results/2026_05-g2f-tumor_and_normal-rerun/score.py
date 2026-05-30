#!/usr/bin/env python3

import argparse
import pandas as pd
from polymerization.score import *
import yaml

parser = argparse.ArgumentParser(description='Calculate fusion scores')
parser.add_argument('--fusion_evidence', default='query-evidence-filled.tsv')
parser.add_argument('--normal_colmap', default='normal_colmap.yaml')
parser.add_argument('--tumor_colmap', default='tumor_colmap.yaml')
parser.add_argument('--output', default='query-scored.tsv')
args = parser.parse_args()

with open(args.normal_colmap) as f:
    normal_colmap = yaml.safe_load(f)
with open(args.tumor_colmap) as f:
    tumor_colmap = yaml.safe_load(f)

df = pd.read_csv(args.fusion_evidence, sep='\t')
low_cov_1000g_cols = ['samples_low_coverage_1000g_dna', 'reads_low_coverage_1000g_dna']
df = df.drop(columns=low_cov_1000g_cols)
# drop low coverage 1000g
cols = df.columns.tolist()
### normal sample columns
print("# scoring normal sample columns")
normal_sample_columns = []
for col in cols:
    if ('1000g' in col) and ('samples' in col):
        normal_sample_columns.append(col)
    if ('normal' in col) and ('samples' in col):
        normal_sample_columns.append(col)
for col in normal_sample_columns:
    for key, sub_dict in normal_colmap.items():
        s_col = sub_dict['samples_col']
        if col == s_col:
            total_samples = sub_dict['total_samples']
            df[col] = sample_score_vectorized(df[col].values, total_samples)
            # rename
            df.rename(columns={col: f'{col}_score'}, inplace=True)
### normal read columns
print("# scoring normal read columns")
normal_read_columns = []
for col in cols:
    if ('1000g' in col) and ('reads' in col):
        normal_read_columns.append(col)
    if ('normal' in col) and ('reads' in col):
        normal_read_columns.append(col)
for col in normal_read_columns:
    for key, sub_dict in normal_colmap.items():
        r_col = sub_dict['reads_col']
        if col == r_col:
            total_samples = sub_dict['total_samples']
            total_upper_bound = sub_dict['upper_bound']
            df[col] = read_score_normal_vectorized(df[col].values, total_samples, total_upper_bound)
            # rename
            df.rename(columns={col: f'{col}_score'}, inplace=True)
### tumor sample columns
print("# scoring tumor sample columns")
tumor_sample_columns = []
for col in cols:
    if ('tumor' in col) and ('samples' in col):
        tumor_sample_columns.append(col)
for col in tumor_sample_columns:
    for key, sub_dict in tumor_colmap.items():
        s_col = sub_dict['samples_col']
        if col == s_col:
            total_samples = sub_dict['total_samples']
            df[col] = sample_score_vectorized(df[col].values, total_samples)
            # rename
            df.rename(columns={col: f'{col}_score'}, inplace=True)
### tumor read columns
print("# scoring tumor read columns")
tumor_read_columns = []
for col in cols:
    if ('tumor' in col) and ('reads' in col):
        tumor_read_columns.append(col)
for col in tumor_read_columns:
    for key, sub_dict in tumor_colmap.items():
        r_col = sub_dict['reads_col']
        if col == r_col:
            total_samples = sub_dict['total_samples']
            total_upper_bound = sub_dict['upper_bound']
            df[col] = read_score_tumor_vectorized(df[col].values, total_samples, total_upper_bound)
            # rename
            df.rename(columns={col: f'{col}_score'}, inplace=True)  

# scale
score_cols = [col for col in df.columns if col.endswith('_score')]
for col in score_cols:
    df[col] = df[col] * 0.5

    
### left off here.
### double check code below

### compute aggregate tumor score
total_tumor_samples = 0
for key, sub_dict in tumor_colmap.items():
    total_tumor_samples += sub_dict['total_samples']
# compute within subpop score
tumor_score_columns = [col for col in df.columns if ('tumor' in col) and col.endswith('_score')]
for tumor_score_col in tumor_score_columns:
    for key, sub_dict in tumor_colmap.items():
        # add read/sample contribution to subpopulation score
        if key in tumor_score_col:
            if key in df.columns:
                df[key] += df[tumor_score_col]
            else:
                df[key] = df[tumor_score_col]
# normalize by total tumor samples
for key, sub_dict in tumor_colmap.items():
    if key in df.columns:
        df[key] = df[key] * (sub_dict['total_samples'] / total_tumor_samples)
### compute aggregate normal score
total_normal_samples = 0
for key, sub_dict in normal_colmap.items():
    total_normal_samples += sub_dict['total_samples']
# compute within subpop score
normal_score_columns = [f'{col}_score' for col in normal_read_columns]
normal_score_columns += [f'{col}_score' for col in normal_sample_columns]
for normal_score_col in normal_score_columns:
    for key, sub_dict in normal_colmap.items():
        # add read/sample contribution to subpopulation score
        if key in normal_score_col:
            if key in df.columns:
                df[key] += df[normal_score_col]
            else:
                df[key] = df[normal_score_col]
# normalize by total normal samples
for key, sub_dict in normal_colmap.items():
    if key in df.columns:
        df[key] = df[key] * (sub_dict['total_samples'] / total_normal_samples)
### compute final tumor and normal scores by summing subpop scores
tumor_subpop_score_columns = [key for key in tumor_colmap.keys() if key in df.columns]
normal_subpop_score_columns = [key for key in normal_colmap.keys() if key in df.columns]
df['tumor_score'] = df[tumor_subpop_score_columns].sum(axis=1)
df['normal_score'] = df[normal_subpop_score_columns].sum(axis=1)
df['final_score'] = df['tumor_score'] - df['normal_score']
df = df.sort_values(by='final_score', ascending=False)
df.to_csv(args.output, sep='\t', index=False)
        

