#!/usr/bin/env python3

import argparse
from polymerization.stix2fusion import left_sort_fusion_set
from polymerization.io import *

parser = argparse.ArgumentParser()
parser.add_argument("--bed", default="../../results/2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed")
parser.add_argument("--fusions",default="fusion_set.tsv")
parser.add_argument("--output",default="fusion_set.sort.tsv")
args = parser.parse_args()

def main():
	df_bed = read_bed(args.bed, gene_col_idx=3)
	df_fusion = read_fusion_set(args.fusions)
	df_fusion_sort = left_sort_fusion_set(df_fusion,df_bed)
	df_fusion_sort.to_csv(args.output,index=False,sep="\t")

if __name__ == "__main__":
    main()



