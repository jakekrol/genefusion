#!/usr/bin/env python3
import argparse
import pandas as pd
import os,sys
from genefusion.genefusion import add_left_right_col
import yaml

parser = argparse.ArgumentParser(description='Setup haas data for 1000 genomes quantification experiment')
parser.add_argument('--haas', type=str,
                    default='../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt',
                    help='Path to haas supplementary table 4')
parser.add_argument('-b','--bedfile', type=str, required=False, help='Path to gene bed file',
                    default='../2025_09-gene_bed/grch37.bed')
parser.add_argument('-a', '--alias_map', type=str, required=False, help='Path to gene alias map',
                    default='../2025_09-haas_interval_frac/merged_alias.yaml')
parser.add_argument('-o','--output', type=str, help='Path to output directory',
                    default='haas_left_right_mapped.tsv')
args = parser.parse_args()

print("Reading input files")
df_haas = pd.read_csv(args.haas, sep='\t', encoding='latin1')
df_bed = pd.read_csv(args.bedfile, sep='\t', header=None,
                     names=['chrom','start','end','gene','strand'])
alias_map = yaml.safe_load(open(args.alias_map))

# build gene set from bed
genes = set(df_bed['gene'].tolist())

# try mapping haas genes to bed gene set using alias map
def map_gene(gene, genes, alias_map):
    if gene in genes:
        return gene
    elif gene in alias_map:
        # will use the first match
        for g in alias_map[gene]:
            if g in genes:
                return g
    return gene

print("Mapping haas genes to bed gene set using alias map")
df_haas['x'] = df_haas['fusion'].apply(lambda x: x.split('--')[0])
df_haas['y'] = df_haas['fusion'].apply(lambda x: x.split('--')[1])
df_haas['x'] = df_haas['x'].apply(lambda x: map_gene(x, genes, alias_map))
df_haas['y'] = df_haas['y'].apply(lambda x: map_gene(x, genes, alias_map))

# make left/right columns
print("Mapping haas genes left/right columns")
df_haas = add_left_right_col(df_haas, 'x', 'y', left_col='left', right_col='right')
cols = df_haas.columns.tolist()
# move left/right to front
cols = ['left','right'] + [c for c in cols if c not in ['left','right']]
df_haas = df_haas[cols]

df_haas.drop(columns=['x','y'], inplace=True)
                    
print("Writing output to", args.output)
df_haas.to_csv(args.output, sep='\t', index=False)



