#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import time
import subprocess

# to-do: verify whether NAs have PE count of 0 or are not part of our gene bed file
LOGREG_EVAL = '/data/jake/genefusion/scripts/python/logreg_eval.py'
parser = argparse.ArgumentParser(description='Run logistic regression evaluation on various random negative fusions')
parser.add_argument('-p', '--positives', type=str, required=True, help='Input positive fusion file')
parser.add_argument('-n', '--negatives', type=str, required=True, help='Input negative fusion file')
parser.add_argument('-d', '--dir-experiment', type=str, required=True, help='Directory for experiment data')
# parser.add_argument('-o', '--output', type=str, required=True, help='Output pattern for results files')
parser.add_argument('-r', '--use_R', action='store_true', help='Use R value in feature set')
parser.add_argument('-m', '--permutations', type=int, default=100, help='Number of permutations to run')
parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed for sampling negatives')
parser.add_argument('-fns', '--false-negatives', type=int, default=None, help='Number of false negatives to use in logreg_eval (default: None)')
parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use for logreg_eval (default: 1)')
args = parser.parse_args()
assert os.path.exists(args.positives), f'Positives file {args.positives} does not exist'
assert os.path.exists(args.negatives), f'Negatives file {args.negatives} does not exist'
assert os.path.exists(args.dir_experiment), f'Experiment directory {args.dir_experiment} does not exist'

# feature cols
stats_cols = ['pe_count', 'sample_count', 'sample_density', 'log10_burden_product']
if args.use_R:
    stats_cols.append('R_transformed')

# read
t = time.time()
df_p = pd.read_csv(args.positives, sep='\t')
print(f'Loaded positives file in {time.time() - t:.2f} seconds')
t = time.time()
df_n = pd.read_csv(args.negatives, sep='\t')
print(f'Loaded negatives file in {time.time() - t:.2f} seconds')
assert all(col in df_p.columns for col in stats_cols), f'Positives file {args.positives} should have columns: {stats_cols}'
assert all(col in df_n.columns for col in stats_cols), f'Negatives file {args.negatives} should have columns: {stats_cols}'

# keep only relevant columns
df_p = df_p[['left', 'right'] + stats_cols]
df_n = df_n[['left', 'right'] + stats_cols]

if args.use_R:
    # some R values were not computed bc of too few data points
    # we replace them with 0 indicating the distance is random
    # apply na_R_to_zero vectorized to R_transformed columns
    df_p['R_transformed'] = df_p['R_transformed'].fillna(0)
    df_n['R_transformed'] = df_n['R_transformed'].fillna(0)

# drop nas in positive fusions
before = df_p.shape[0]
df_p = df_p.dropna()
print(f'Dropped {before - df_p.shape[0]} rows with NAs in positives')

# drop nas in negative fusions
before = df_n.shape[0]
df_n = df_n.dropna()
print(f'Dropped {before - df_n.shape[0]} rows with NAs in negatives')
assert before == df_n.shape[0], 'Negatives should not have NAs, but they do. Inspect w/ breakpoint()'

# add label if not present
if 'label' not in df_p.columns:
    df_p['label'] = 1  # positive samples
if 'label' not in df_n.columns:
    df_n['label'] = 0  # negative samples

# make args.permutation dataset mixtures
for i in range(args.permutations):
    print(f'Processing permutation {i}/{args.permutations-1}')
    # sample negatives
    df_n_rand = df_n.sample(n=df_p.shape[0], random_state=args.seed + i) # add i to seed to ensure different, yet consistent samples each time
    # combine positives and negatives
    df_combined = pd.concat([df_p, df_n_rand], ignore_index=True)
    # write the fusions
    fusions_out = os.path.join(args.dir_experiment, f'rand_{i}_fusions.tsv')
    df_combined[['left', 'right']].to_csv(fusions_out, sep='\t', index=False)
    # write the feature and label data
    rand_out = os.path.join(args.dir_experiment, f'rand_{i}_features.tsv')
    df_combined.drop(columns=['left', 'right'], inplace=True)  # drop fusion columns
    df_combined.to_csv(rand_out, sep='\t', index=False)
    # get label idx
    if 'label' in df_combined.columns:
        label_idx = df_combined.columns.get_loc('label') + 1  # +1 for 1-indexed in logreg_eval
    else:
        raise ValueError('Label column not found in combined DataFrame')
    # run logreg_eval
    cmd = [
        LOGREG_EVAL,
        '-i', rand_out,
        '-y', str(label_idx),
        '-k', '5',  # 5-fold CV
        '-o', os.path.join(args.dir_experiment, f'rand_{i}_logreg'),
        '--seed', str(args.seed + i),
        '--cpus', str(args.cpus),
    ]
    if args.false_negatives:
        cmd += ['-fns', str(args.false_negatives)]
    t = time.time()
    print(f'Running logreg_eval for permutation {i} with command: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    print(f'Ran logreg_eval for permutation {i} in {time.time() - t:.2f} seconds')


    

