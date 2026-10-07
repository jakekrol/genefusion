#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

parser=argparse.ArgumentParser()
parser.add_argument("--input",default="fusion_eval.scored.tsv")
parser.add_argument("--output_datset_overlap", default="dataset_overlap.png")
parser.add_argument("--output_score_hist",default="score-hist.png")
args = parser.parse_args()

score_coverage='score_global_norm_coverage'
score_burden='score_global_norm_burden'

def main():
    df = pd.read_csv(args.input,sep="\t",index_col=0)
    dataset2fusions=defaultdict(list)
    fusion2score_burden=dict()
    fusion2score_coverage=dict()
    fusions = set()
    for i,row in df.iterrows():
        fusion = row['gene_left'] + '--' + row['gene_right']
        fusions.add(fusion)
        fusion2score_burden[fusion] = row[score_burden]
        fusion2score_coverage[fusion] = row[score_coverage]
        for dataset in row['eval_dataset'].split(','):
            dataset2fusions[dataset].append(fusion)
    print(f"# total unique fusions: {len(fusions)}")
    # hist score per dataset
    datasets=dataset2fusions.keys()
    dataset2scores_burden=defaultdict(list)
    for dataset,fusions in dataset2fusions.items():
        for fusion in fusions:
            score = fusion2score_burden[fusion]
            dataset2scores_burden[dataset].append(score)
    dataset2scores_coverage=defaultdict(list)
    for dataset,fusions in dataset2fusions.items():
        for fusion in fusions:
            score = fusion2score_coverage[fusion]
            dataset2scores_coverage[dataset].append(score)

    fig, ax = plt.subplots(2, 6, figsize=(18, 6))

    for i, (k, v) in enumerate(dataset2scores_burden.items()):
        ax[0, i].hist(v, color='black')
        ax[0, i].set_title(k, fontsize=9, loc='left')
        ax[0, i].set_yscale('log')
        ax[0, i].set_xlim(-1, 1)
        ax[0, i].spines['top'].set_visible(False)
        ax[0, i].spines['right'].set_visible(False)

    for i, (k, v) in enumerate(dataset2scores_coverage.items()):
        ax[1, i].hist(v, color='black')
        ax[1, i].set_title(k, fontsize=9, loc='left')
        ax[1, i].set_yscale('log')
        ax[1, i].set_xlim(-1, 1)
        ax[1, i].spines['top'].set_visible(False)
        ax[1, i].spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(args.output_score_hist)
    plt.close()

    datasets = list(dataset2fusions)
    overlap = []

    for d1 in datasets:
        row = []
        f1 = set(dataset2fusions[d1])
        for d2 in datasets:
            row.append(len(f1 & set(dataset2fusions[d2])))
        overlap.append(row)

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(overlap, cmap='Reds')

    ax.set_xticks(range(len(datasets)), datasets, rotation=90)
    ax.set_yticks(range(len(datasets)), datasets)

    ax.set_title('Fusion overlap', loc='left', fontsize=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for i in range(len(datasets)):
        for j in range(len(datasets)):
            ax.text(j, i, overlap[i][j], ha='center', va='center', fontsize=8)

    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(args.output_datset_overlap)
    plt.close()
    

if __name__ == "__main__":
    main()