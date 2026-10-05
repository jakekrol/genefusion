#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
from polymerization.analysis import intersect_file2sample_counts
# from polymerization.datasets import get_pcawg_data_types

parser=argparse.ArgumentParser()
parser.add_argument('--indir',default='top_fusions')
parser.add_argument('--outdir',default='top_fusions')
parser.add_argument('--infile_suffix', default=".annotate.repeat_filt.breakpoint.rank.tsv")
parser.add_argument("--outfile_suffix", default= ".annotate.repeat_filt.breakpoint.rank.samples.tsv")
parser.add_argument("--num_fusions", default=10, help="number of fusions to analyze per file")
parser.add_argument("--num_samples", default=3, help="number of samples with maximal read evidence to analyze")
parser.add_argument("--dir_g2f",default="../2026_06-g2f-all_gene_pairs/g2f_out")
parser.add_argument("--intersect_file_suffix",default=".giggle.clean.swap.intersect.bed.gz")
args = parser.parse_args()
args.num_fusions=int(args.num_fusions)
args.num_samples=int(args.num_samples)

# focus on tissue with rna-seq
TISSUES={'blood', 'kidney', 'liver', 'ovary'}
MODALITY='rna'

def main():
    files = glob.glob(f"{args.indir}/*{args.infile_suffix}")
    for f in files:
        for tissue in TISSUES:
            if tissue in f:
                # find g2f subdir
                for directory in os.listdir(args.dir_g2f):
                    if (tissue in directory) and (MODALITY in directory) and ("tumor" in directory):
                        g2f_subdir = directory
                df = pd.read_csv(f,sep="\t")
                df = df.head(args.num_fusions)
                df['top_samples'] = pd.NA
                j=0
                for i,row in df.iterrows():
                    gene_left = row['gene_left']
                    gene_right = row['gene_right']
                    intersect_file = os.path.join(
                        os.path.join(
                            args.dir_g2f,
                            g2f_subdir,
                            f"{gene_left}{args.intersect_file_suffix}"
                        )
                    )
                    if os.path.exists(intersect_file):
                        sample_counts = intersect_file2sample_counts(intersect_file, gene_right, bgzip=True)
                        top_samples = sorted(sample_counts, key=sample_counts.get, reverse=True)[:args.num_samples]
                        # extract pcawg file id from sample name string
                        # format samples as csv string
                        top_samples = ','.join([os.path.basename(x).split('.')[0] for x in top_samples])
                        df.at[i,'top_samples'] = top_samples
                outfile = f.replace(args.infile_suffix,args.outfile_suffix)
                df.to_csv(outfile,sep="\t",index=False)

if __name__ == "__main__":
    main()