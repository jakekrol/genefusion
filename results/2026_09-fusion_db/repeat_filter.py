#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--repeat_coverage",default = "../2026_09-gene_repeat_coverage/gene_repeat.coverage.tsv")
parser.add_argument("--indir",default="top_fusions")
parser.add_argument("--infile_suffix", default=".annotated.tsv")
parser.add_argument("--outfile_suffix", default=".annotate.repeat_filt.tsv")
parser.add_argument("--outdir", default="top_fusions")
parser.add_argument("--threshold", default=0.75, help="repeat coverage upper bound inclusion")
args = parser.parse_args()
args.threshold=float(args.threshold)

def main():
    df_cov = pd.read_csv(args.repeat_coverage,sep="\t")
    gene2rcov=dict()
    for i,row in df_cov.iterrows():
        gene = row['gene_name']
        repeat_cov = row['fraction_of_gene_intersecting_repeat']
        gene2rcov[gene] = repeat_cov
    files = glob.glob(f"{args.indir}/*{args.infile_suffix}")
    for f in files:
        df_fusion = pd.read_csv(f,sep='\t')
        df_fusion.columns=['fusion','score','annotation']
        df_fusion['gene_left'] = df_fusion['fusion'].apply(lambda x: x.split('--')[0])
        df_fusion['gene_right'] = df_fusion['fusion'].apply(lambda x: x.split('--')[1])
        df_fusion['repeat_cov_left'] = df_fusion['gene_left'].map(gene2rcov)
        df_fusion['repeat_cov_right'] = df_fusion['gene_right'].map(gene2rcov)
        mask = (df_fusion['repeat_cov_left'] <= args.threshold) & \
            (df_fusion['repeat_cov_right'] <= args.threshold)
        df_fusion=df_fusion[mask].reset_index(drop=True)
        out = os.path.join(
            args.outdir,
            os.path.basename(f).replace(args.infile_suffix, args.outfile_suffix)
        )
        df_fusion.to_csv(out,sep="\t",index=False)
    
if __name__ == "__main__":
    main()