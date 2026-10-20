#!/usr/bin/env python3
import argparse
import ast
import os
import pandas as pd
import subprocess
from polymerization.analysis import intersect_file2breakpoints, estimate_breakpoint
from polymerization.score import score_breakpoint
GENE_LEFT="RP11-345J4.5"
GENE_RIGHT="RP11-345J4.8"
parser = argparse.ArgumentParser()
parser.add_argument("--samples", default="samples.tsv")
parser.add_argument("--exons", default= "../2026_09-grch37-gtf2exon_intervals/exon_intervals.tsv.gz")
parser.add_argument("--default_score",default=2**32)
parser.add_argument(
    "--intersect",
    default=f"../2026_06-g2f-all_gene_pairs/g2f_out/kidney_normal_rna/{GENE_LEFT}.giggle.clean.swap.intersect.bed.gz"
)
parser.add_argument("--output", default="breakpoints.tsv")
args = parser.parse_args()

def get_subset_intersect(file_intersect: str, samples: list):
    df = intersect_file2breakpoints(file_intersect, GENE_RIGHT, True, True)
    mask = df['sample'].isin(samples)
    df = df[mask].reset_index(drop=True)
    return df


def main():
    df_exon = pd.read_csv(args.exons,sep="\t")
    mask = df_exon['gene'].isin([GENE_LEFT,GENE_RIGHT])
    df_exon = df_exon[mask].reset_index(drop=True)
    df_exon.set_index('gene',inplace=True)
    df_samples = pd.read_csv(args.samples,sep="\t")
    samples = df_samples['sample'].tolist()
    df_intersect_subset = get_subset_intersect(args.intersect, samples)
    bp_left, bp_right = estimate_breakpoint(df_intersect_subset)
    df_samples['breakpoint_left'] = bp_left
    df_samples['breakpoint_right'] = bp_right
    transcript_left,_ = score_breakpoint(GENE_LEFT, df_exon, int(bp_left.split(':')[1]),args.default_score)
    transcript_right,_ = score_breakpoint(GENE_RIGHT, df_exon, int(bp_right.split(':')[1]),args.default_score)
    df_samples['transcript_left'] = transcript_left
    df_samples['transcript_right'] = transcript_right
    df_samples.to_csv(args.output,sep="\t",index=False)

if __name__ == "__main__":
    main()
