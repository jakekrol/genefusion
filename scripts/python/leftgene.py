#!/usr/bin/env python3
import pandas as pd
import os,sys
import argparse

# input: 2 genes, and gene bed file
# output: gene with smaller position on the genome

parser = argparse.ArgumentParser(description='Get the gene with smaller position on the genome')
parser.add_argument('-x', '--genex', type=str, required=True, help='Gene x')
parser.add_argument('-y', '--geney', type=str, required=True, help='Gene y')
parser.add_argument('-b', '--bedfile', type=str, default='/data/jake/genefusion/data/gene_file.txt.latest')
args = parser.parse_args()

df = pd.read_csv(args.bedfile, sep='\t', header=None)

def get_gene_position(gene, df):
    # Filter the dataframe for the gene
    df_gene = df[df[3] == gene]
    if df_gene.shape[0] == 0:
        print(-1)
        sys.exit(1)
    # Get the position of the gene
    # bug: need to handle X and Y
    # chrom = int(df_gene.iloc[0, 0]) 
    chrom = df_gene.iloc[0, 0]
    start = int(df_gene.iloc[0, 1])
    end = int(df_gene.iloc[0, 2])
    return chrom, start, end
def compare_genes(gene1, gene2, df):
    chrom1, start1, end1 = get_gene_position(gene1, df)
    chrom2, start2, end2 = get_gene_position(gene2, df)
    # handle X and Y chromosomes
    # both are sex chromosomes
    if (chrom1 in ['X', 'Y']) and (chrom2 in ['X', 'Y']):
        # compare X and Y chromosomes
        if chrom1 == 'X' and chrom2 == 'Y':
            return gene1
        elif chrom1 == 'Y' and chrom2 == 'X':
            return gene2
    # one is sex chromosome and the other is not
    if (chrom1 in ['X', 'Y']) and (chrom2 not in ['X', 'Y']):
        return gene2
    if (chrom2 in ['X', 'Y']) and (chrom1 not in ['X', 'Y']):
        return gene1
    # same sex chromosome
    if (chrom1 == chrom2) and (chrom1 in ['X', 'Y']):
         # start
        if start1 < start2:
            return gene1
        elif start1 > start2:
            return gene2
        else:
            # end
            if end1 < end2:
                return gene1
            else:
                return gene2
    # autosomes
    if chrom1 < chrom2:
        return gene1
    elif chrom1 > chrom2:
        return gene2
    else:
        # start
        if start1 < start2:
            return gene1
        elif start1 > start2:
            return gene2
        else:
            # end
            if end1 < end2:
                return gene1
            else:
                return gene2
left = compare_genes(args.genex, args.geney, df)
print(left.strip())


