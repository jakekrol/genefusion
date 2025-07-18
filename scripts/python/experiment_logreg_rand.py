#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import time
import subprocess

LOGREG_EVAL = '/data/jake/rl-tools/ml/logreg_eval.py'
# args
parser = argparse.ArgumentParser(description='Run logistic regression evaluation on various random negative fusions')
parser.add_argument('-p', '--positives', type=str, required=True, help='Input positive fusion file')
parser.add_argument('-r', '--randoms', type=str, required=True, help='Input random negative fusion file')
parser.add_argument('-d', '--dir-experiment', type=str, required=True, help='Directory for experiment data')
parser.add_argument('-n', '--negatives', type=str, default=None, help='Input negative fusion file')
parser.add_argument('-m', '--permutations', type=int, default=100, help='Number of permutations to run')
parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed for sampling negatives')
parser.add_argument('-fns', '--false-negatives', type=int, default=None, help='Number of false negatives to use in logreg_eval (default: None)')
parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use for logreg_eval (default: 1)')
parser.add_argument('--tune_metric', type=str, default='auroc', choices=['auroc', 'auprc'], help='Metric to tune for logreg_eval (default: auroc)')
args = parser.parse_args()
if args.negatives is None:
    print('No negatives file specified, using randoms as negatives:', args.randoms)
assert os.path.exists(args.positives), f'Positives file {args.positives} does not exist'
assert os.path.exists(args.randoms), f'Randoms file {args.randoms} does not exist'
assert os.path.exists(args.dir_experiment), f'Experiment directory {args.dir_experiment} does not exist'

# read
t = time.time()
df_p = pd.read_csv(args.positives, sep='\t')
print(f'Loaded positives file in {time.time() - t:.2f} seconds')
t = time.time()
df_r = pd.read_csv(args.randoms, sep='\t')
print(f'Loaded randoms file in {time.time() - t:.2f} seconds')
if args.negatives:
    t = time.time()
    df_n = pd.read_csv(args.negatives, sep='\t')
    print(f'Loaded negatives file in {time.time() - t:.2f} seconds')

# drop nas in positive fusions
before = df_p.shape[0]
df_p = df_p.dropna()
print(f'Dropped {before - df_p.shape[0]} rows with NAs in positives')

# drop nas in random fusions
before = df_r.shape[0]
df_r = df_r.dropna()
print(f'Dropped {before - df_r.shape[0]} rows with NAs in randoms')
assert before == df_r.shape[0], 'Randoms should not have NAs, but they do. Inspect w/ breakpoint()'

# drop nas in negatives if provided
if args.negatives:
    before = df_n.shape[0]
    df_n = df_n.dropna()
    print(f'Dropped {before - df_n.shape[0]} rows with NAs in negatives')
    assert before == df_n.shape[0], 'Negatives should not have NAs, but they do. Inspect w/ breakpoint()'

# add label if not present
if 'label' not in df_p.columns:
    df_p['label'] = 1  # positive samples
if 'label' not in df_r.columns:
    df_r['label'] = 0  # negative samples
if args.negatives and 'label' not in df_n.columns:
    df_n['label'] = 0  # negative samples

print(f'Positives shape: {df_p.shape}'
      f'\nRandoms shape: {df_r.shape}'
      f'\nNegatives shape: {df_n.shape if args.negatives else "N/A"}')
if args.negatives:
    n_samples = df_p.shape[0] - df_n.shape[0]
else:
    n_samples = df_p.shape[0]
print(f'Number of samples to draw from randoms: {n_samples}')

# make args.permutation dataset mixtures
for i in range(args.permutations):
    print(f'Processing permutation {i}/{args.permutations-1}')
    # sample negatives
    df_i = df_r.sample(n=n_samples, random_state=args.seed + i) # add i to seed to ensure different, yet consistent samples each time
    # stack
    if args.negatives:
        df_combined = pd.concat([df_p, df_i, df_n], ignore_index=True)
    else:
        df_combined = pd.concat([df_p, df_i], ignore_index=True)
    # write fusions and feature data
    i_out_fusion = os.path.join(args.dir_experiment, f'rand_{i}_fusions.tsv')
    df_combined[['left', 'right', 'label']].to_csv(i_out_fusion, sep='\t', index=False)
    df_combined = df_combined.drop(columns=['left', 'right'])  # drop fusion columns to keep only features and label
    i_out_data = os.path.join(args.dir_experiment, f'rand_{i}_features.tsv')
    df_combined.to_csv(i_out_data, sep='\t', index=False)
    # get label idx
    if 'label' in df_combined.columns:
        label_idx = df_combined.columns.get_loc('label') + 1  # +1 for 1-indexed in logreg_eval
    else:
        raise ValueError('Label column not found in combined DataFrame')
    # run logreg_eval
    cmd = [
        LOGREG_EVAL,
        '-i', i_out_data,
        '-y', str(label_idx),
        '-k', '5',  # 5-fold CV
        '-o', os.path.join(args.dir_experiment, f'rand_{i}_logreg'),
        '--seed', str(args.seed + i),
        '--cpus', str(args.cpus),
        '--tune_metric', args.tune_metric
    ]
    if args.false_negatives:
        cmd += ['-fns', str(args.false_negatives)]
    t = time.time()
    print(f'Running logreg_eval for permutation {i} with command: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    print(f'Ran logreg_eval for permutation {i} in {time.time() - t:.2f} seconds')


    

