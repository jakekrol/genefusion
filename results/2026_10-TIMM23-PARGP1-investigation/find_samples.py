#!/usr/bin/env python3
import argparse
import pandas as pd
from polymerization.analysis import intersect_file2breakpoints
parser = argparse.ArgumentParser()
parser.add_argument("--thousg_rna",default="../2026_06-g2f-all_gene_pairs/g2f_out/mage_short_read_1000g_rna/PARGP1.giggle.clean.swap.intersect.bed.gz")
parser.add_argument("--pcawg_blood_tumor",default="../2026_06-g2f-all_gene_pairs/g2f_out/blood_tumor_rna/PARGP1.giggle.clean.swap.intersect.bed.gz")
parser.add_argument("--output",default="samples.tsv")
args = parser.parse_args()

GENE_RIGHT='TIMM23'
COL_SAMPLE='sample'
K=3

def main():
	df_thousg_rna = intersect_file2breakpoints(args.thousg_rna, GENE_RIGHT, False, True)
	s_thousg_rna = df_thousg_rna[COL_SAMPLE].value_counts().head(K)
	df_pcawg_blood_tumor = intersect_file2breakpoints(args.pcawg_blood_tumor, GENE_RIGHT, False, True)
	s_pcawg_blood_tumor = df_pcawg_blood_tumor[COL_SAMPLE].value_counts().head(K)
	ser = pd.concat([s_thousg_rna, s_pcawg_blood_tumor])
	ser.to_csv(args.output,sep="\t")

if __name__ == "__main__":
    main()

