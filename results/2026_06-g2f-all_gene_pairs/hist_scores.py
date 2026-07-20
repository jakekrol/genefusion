#!/usr/bin/env python3

import argparse
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from polymerization.datasets import get_pcawg_recurrent_tumor_fusions
from polymerization.datasets import get_recurrent_normal_tissue_specific_fusions

parser = argparse.ArgumentParser(description='Tissue-wise score hists')
parser.add_argument('--score', required=True)
parser.add_argument('--seed', type=int, default=0)
parser.add_argument('--sample_size', type=int, default=100000)
parser.add_argument('--out_data', required=True)
parser.add_argument('--out_png', required=True)
parser.add_argument('--bins', type=int, default=50)
parser.add_argument('--title', default='')
parser.add_argument('--tissue', required=True)
args = parser.parse_args()

def lookup_fusion_score(gene_left, gene_right, score_file):
    cmd = ["grep", rf"^{gene_left}\s{gene_right}\s", score_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        score = result.stdout.split('\t')[2].strip()
        return float(score)
    else:
        return None

print("# gathering recurrent tumor fusions for tissue", args.tissue)
# tumor fusions
df_tumor = get_pcawg_recurrent_tumor_fusions()
mask = df_tumor['tissues'].apply(lambda x: args.tissue.lower() in x.lower().split(','))
df_tumor = df_tumor[mask]
tumor_fusions = []
for i,row in df_tumor.iterrows():
    gene_left = row['gene_left']
    gene_right = row['gene_right']
    score = None
    score = lookup_fusion_score(gene_left, gene_right, args.score)
    if score is not None:
        tumor_fusions.append((gene_left, gene_right, score))
# normal fusions
print("# gathering recurrent normal fusions for tissue", args.tissue)
df_normal = get_recurrent_normal_tissue_specific_fusions()
mask = df_normal['tissues'].apply(lambda x: args.tissue.lower() in x.lower().split(','))
df_normal = df_normal[mask]
normal_fusions = []
for i,row in df_normal.iterrows():
    gene_left = row['gene_left']
    gene_right = row['gene_right']
    score = None
    score = lookup_fusion_score(gene_left, gene_right, args.score)
    if score is not None:
        normal_fusions.append((gene_left, gene_right, score))

    
    

print("# sampling with shuf using seed", args.seed)

# Use shuf with openssl deterministic random source
# cmd = (
#     "tail -n +2 {score} | "
#     "shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:{seed} -nosalt </dev/zero 2>/dev/null) "
#     "-n {sample_size} | "
#     "cut -f3 > {out_data}"
# ).format(
#     score=args.score,
#     seed=args.seed,
#     sample_size=args.sample_size,
#     out_data=args.out_data
# )

# subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

print("# plotting histogram")
x = []
with open(args.out_data, 'r') as f:
	for line in f:
		x.append(float(line.strip()))

# Split into 3 regions
region1 = [v for v in x if -1 <= v < -0.1]  # [-1, -0.1)
region2 = [v for v in x if -0.1 <= v < 0.1]  # [-0.1, 0.1)
region3 = [v for v in x if 0.1 <= v <= 1]   # [0.1, 1]

# find which region the tumor and normal fusions fall into
r1_tumor = []
r2_tumor = []
r3_tumor = []
r1_normal = []
r2_normal = []
r3_normal = []
for gene_left, gene_right, score in tumor_fusions:
    if -1 <= score < -0.1:
        r1_tumor.append(score)
    elif -0.1 <= score < 0.1:
        r2_tumor.append(score)
    elif 0.1 <= score <= 1:
        r3_tumor.append(score)
    else:
        pass
for gene_left, gene_right, score in normal_fusions:
    if -1 <= score < -0.1:
        r1_normal.append(score)
    elif -0.1 <= score < 0.1:
        r2_normal.append(score)
    elif 0.1 <= score <= 1:
        r3_normal.append(score)
    else:
        pass

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
fig.suptitle(args.title, fontsize=14)

# Region 1: [-1, -0.1)
axes[0].hist(region1, bins=args.bins, color='black', edgecolor='black')
axes[0].set_xlabel('Score')
axes[0].set_ylabel('Count')
axes[0].set_title(f'[-1, -0.1) (n={len(region1)})')
axes[0].set_xlim([-1, -0.1])
for score in r1_tumor:
    axes[0].axvline(score, color='orange', linestyle='--', linewidth=1)
for score in r1_normal:
    axes[0].axvline(score, color='blue', linestyle='--', linewidth=1)

# Region 2: [-0.1, 0.1)
axes[1].hist(region2, bins=args.bins, color='black', edgecolor='black')
axes[1].set_yscale('log')
axes[1].set_xlabel('Score')
axes[1].set_ylabel('log10(Count)')
axes[1].set_title(f'[-0.1, 0.1) (n={len(region2)})')
axes[1].set_xlim([-0.1, 0.1])
for score in r2_tumor:
    axes[1].axvline(score, color='orange', linestyle='--', linewidth=1)
for score in r2_normal:
    axes[1].axvline(score, color='blue', linestyle='--', linewidth=1)

# Region 3: [0.1, 1]
axes[2].hist(region3, bins=args.bins, color='black', edgecolor='black')
axes[2].set_xlabel('Score')
axes[2].set_ylabel('Count')
axes[2].set_title(f'[0.1, 1] (n={len(region3)})')
axes[2].set_xlim([0.1, 1])
for score in r3_tumor:
    axes[2].axvline(score, color='orange', linestyle='--', linewidth=1)
for score in r3_normal:
    axes[2].axvline(score, color='blue', linestyle='--', linewidth=1)

plt.tight_layout()
plt.savefig(args.out_png, dpi=150)
print(f"# saved plot to {args.out_png}")
