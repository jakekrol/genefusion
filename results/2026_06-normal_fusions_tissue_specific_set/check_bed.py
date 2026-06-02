#!/usr/bin/env python3

from polymerization.io import *
import argparse

# check if bedfile has the recurrent normal tissue specific fusions

parser = argparse.ArgumentParser(description='Check if bedfile has the recurrent normal tissue specific fusions')
parser.add_argument('--fusions', default = './recurrent_normal_fusion_set.tsv')
parser.add_argument('--bedfile', default = '../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed', help='input bed file')
parser.add_argument('--outfile', default = 'missing_bed_genes.txt', help='output file')
args = parser.parse_args()

df_bed = read_bed(args.bedfile, gene_col_idx=3)
df_fusion  = pd.read_csv(args.fusions, sep='\t')
A = set(df_fusion['up_gene'].unique())
B = set(df_fusion['dw_gene'].unique())
fusion_genes = A.union(B)
genes_in_bed = set(df_bed['gene_name'].unique())
genes_not_in_bed = fusion_genes - genes_in_bed
n = len(fusion_genes)
m = len(genes_not_in_bed)
with open(args.outfile, 'w') as f:
    for gene in genes_not_in_bed:
        f.write(gene + '\n')
print(f"{m} out of {n} fusion genes are missing from the bed file. See '{args.outfile}' for details.")


