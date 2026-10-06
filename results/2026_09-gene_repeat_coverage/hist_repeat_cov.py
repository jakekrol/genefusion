#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("--repeat_bed_coverage",default="gene_repeat.coverage.bed")
parser.add_argument("--output", default="gene-repeat-coverage.hist.png")
parser.add_argument("--title", default="GRCh37 genes")
parser.add_argument("--bins", default=30)
parser.add_argument("--threshold",default=0.75)
args = parser.parse_args()
args.threshold=float(args.threshold)

def plot(x,text):
    fig, ax = plt.subplots(figsize=(6,4))
    ax.hist(x, color='black',bins=args.bins)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(args.title,loc='left')
    ax.set_xlabel("Gene repeat coverage")
    ax.set_ylabel("Count")
    ax.axvline(args.threshold, linestyle='--', color='red')
    ax.annotate(text,xy=(0.15,0.8),xycoords="axes fraction")
    fig.savefig(args.output)

def main():
    df = pd.read_csv(
        args.repeat_bed_coverage,
        sep="\t"
    )
    df.columns = [
        'chrom', 'start','end', 'gene_name', 'strand',
        'num_repeat_intervals_intersecting_gene', 
        'num_bases_gene_intersecting_repeat',
        'gene_length',
        'fraction_of_gene_intersecting_repeat'
    ]
    m = df.shape[0]
    mask = df['fraction_of_gene_intersecting_repeat'] > args.threshold
    frac_gt_half = round(mask.sum()/m, 3)
    print(f"# {frac_gt_half} genes have >{args.threshold} repeat coverage")
    text = f"{frac_gt_half* 100}% genes > {args.threshold} repeat coverage"
    plot(df['fraction_of_gene_intersecting_repeat'],text)
    out_data = args.repeat_bed_coverage.replace(".bed", ".tsv")
    df.to_csv(out_data,sep="\t", index=False)

    
if __name__ == "__main__":
    main()