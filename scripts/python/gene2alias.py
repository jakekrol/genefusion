#!/usr/bin/env python3
import os,sys
import pandas as pd
import argparse
import yaml
import swifter

# map source genes to targets using alias file
parser = argparse.ArgumentParser(description='Map gene names to aliases')
parser.add_argument('--source', type=str, required=True, help='Input file with gene names')
parser.add_argument('--source_cols', type=str, help='Comma seprated column indices (1-based) in source file to map', default='1')
parser.add_argument('--source_header', action='store_true', help='Whether source file has a header row')
parser.add_argument('--alias', type=str, required=True, help='Gene alias mapping file') 
parser.add_argument('--target', type=str, required=True, help='Line delimited file of target gene names')
parser.add_argument('--output', type=str, required=True, help='Output file with mapped gene names')
args = parser.parse_args()

# get source columns to map
source_cols = [int(x)-1 for x in args.source_cols.split(',')]

print('loading input files...', file=sys.stderr)
# load source file
source_df = pd.read_csv(args.source, sep='\t', header=0 if args.source_header else None, dtype=str, usecols=source_cols)

# load alias file
with open(args.alias) as f:
    alias_map = yaml.safe_load(f)

# load target genes
target_genes = set()
with open(args.target) as f:
    for line in f:
        target_genes.add(line.strip().upper())

# map source genes to target genes using alias
def map_gene(gene):
    gene_upper = gene.upper()
    if gene_upper in target_genes:
        return gene
    if gene_upper in alias_map:
        for alias in alias_map[gene_upper]:
            if alias.upper() in target_genes:
                return alias
    return gene  # return original if no mapping found

# apply mapping to source dataframe
print('mapping genes...', file=sys.stderr)
for col in source_cols:
    source_df.iloc[:, col] = source_df.iloc[:, col].swifter.apply(map_gene)

# save output
source_df.to_csv(args.output, sep='\t', index=False)
