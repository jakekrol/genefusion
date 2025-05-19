#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

# details:
# visualize the distribution of PE stats for true positives and random negatives

parser = argparse.ArgumentParser(description='Plot histograms of PE stats')
parser.add_argument('-p', '--true_positives', type=str, required=True, help='path to true positives stats table')
parser.add_argument('-n', '--random_negatives', type=str, required=True, help='path to random negatives stats table')
parser.add_argument('-o', '--output', type=str, required=True, help='Output pattern for histograms')
parser.add_argument('-s', '--seed', type=int, default=42, help='Random seed for reproducibility')
args = parser.parse_args()

df_p = pd.read_csv(args.true_positives, sep='\t')
df_n = pd.read_csv(args.random_negatives, sep='\t')
print(f'True Positives shape: {df_p.shape}')
print(f'Random Negatives shape: {df_n.shape}')

# drop sort_result -1 rows
df_p = df_p[df_p['sort_result'] != -1]
print(f'True Positives shape after dropping sort_result=-1: {df_p.shape}')

# add burden product column
df_p['burden_product'] = df_p['burden_total_left'] * df_p['burden_total_right']
df_n['burden_product'] = df_n['burden_total_left'] * df_n['burden_total_right']
df_p['burden_product_log10'] = df_p['burden_product'].apply(lambda x: 0 if x == 0 else np.log10(x))
df_n['burden_product_log10'] = df_n['burden_product'].apply(lambda x: 0 if x == 0 else np.log10(x))

# add sample density column
# in [0,1]
# -> 1 means PE evidence ~ # of samples
# -> 0 means PE evidence >> # of samples
def sample_density(row):
    if row['pe_count'] == 0:
        return -1
    elif row['pe_count'] < row['sample_count']:
        return -1
    else:
        return row['sample_count'] / row['pe_count']
df_p['sample_density'] = df_p.apply(lambda row: sample_density(row), axis=1)
df_n['sample_density'] = df_p.apply(lambda row: sample_density(row), axis=1)

plot_cols = ['pe_count', 'sample_count', 'burden_total_left', 'burden_total_right', 'burden_product_log10', 'sample_density']

# downsample the random negatives to match the number of true positives
m = df_p.shape[0]
df_n = df_n.sample(n=m, random_state=args.seed)

# do an overlay histogram of the two distributions
# for each column in plot_cols
print(f'Plotting histograms for columns: {plot_cols}')
# print the shape of the dataframes
print(f'True Positives shape: {df_p.shape}')
print(f'Random Negatives shape: {df_n.shape}')
for col in plot_cols:
    plt.figure(figsize=(10, 6))
    plt.hist(df_p[col], bins=50, alpha=0.5, label='True Positives', color='blue')
    plt.hist(df_n[col], bins=50, alpha=0.5, label='Random Negatives', color='red')
    plt.xlabel(col)
    plt.ylabel('Frequency')
    plt.title(f'Histogram of {col}')
    plt.legend()
    plt.savefig(f'{args.output}_{col}.png')
    plt.close()


