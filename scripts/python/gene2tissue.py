#!/usr/bin/env python3

import os,sys
import pandas as pd
import argparse
import yaml
# import swifter

parser = argparse.ArgumentParser(description="Map genes to tissues using Human Protein Atlas data")
parser.add_argument("-i", "--input", required=True, help="Input file with gene list (one gene per line)")
parser.add_argument("-d", "--data",
                    default="/data/jake/genefusion/data/2025_10-human_protein_atlas_gene_rna_expr/rna_tissue_consensus.tsv",
                    help="Path to rna_tissue_consensus.tsv file from Human Protein Atlas")
parser.add_argument("-a", "--alias", help="Optional gene alias mapping file",
                    default="/data/jake/genefusion/results/2025_09-haas_interval_frac/merged_alias.yaml")
parser.add_argument("-o", "--out", required=True, help="Output file for gene-tissue mapping")
parser.add_argument("-k", '--topk', type=int, default=3, help="Number of top tissue(s) to report per gene")
parser.add_argument("--out_header", type =str, default=None, help="Custom header for output file (comma-separated)")
# parser.add_argument("--tpm", action="store_true", help="Report TPM values along with tissue(s)")

args = parser.parse_args()

COLUMN_GENE = 'Gene name'

# load gene list
print("Loading gene list from", args.input)
queries = set()
with open(args.input) as f:
    for line in f:
        queries.add(line.strip().upper())
    
# load alias mapping if provided
print("Loading alias mapping from", args.alias)
alias_map = {}
if args.alias:
    with open(args.alias) as f:
        alias_map = yaml.safe_load(f)
    
# load HPA data
print("Loading HPA data from", args.data)
hpa = pd.read_csv(args.data, sep="\t")

def map_gene_to_tissues(gene):
    x = None
    gene_upper = gene.upper()
    # check if gene in hpa
    if gene_upper in hpa[COLUMN_GENE].values:
        mask = hpa[COLUMN_GENE] == gene_upper
        x = hpa[mask][['Tissue', 'nTPM']].sort_values(by='nTPM', ascending=False)
        x = x.head(args.topk)['Tissue']
    # otherwise check alias map
    else:
        if gene_upper in alias_map:
            for alias in alias_map[gene_upper]:
                alias_upper = alias.upper()
                if alias_upper in hpa[COLUMN_GENE].values:
                    mask = hpa[COLUMN_GENE] == alias_upper
                    x = hpa[mask][['Tissue', 'nTPM']].sort_values(by='nTPM', ascending=False)
                    x = x.head(args.topk)['Tissue']
                    break
    if x is None or x.empty:
        return("NA")
    else:
        return(",".join(x.tolist()))

df_q = pd.DataFrame(list(queries), columns=['gene'])
# df_q['Tissues'] = df_q['gene'].swifter.apply(map_gene_to_tissues)
df_q['Tissues'] = df_q['gene'].apply(map_gene_to_tissues)
df_q = df_q.sort_values(by='gene')
if args.out_header:
    df_q.columns = [args.out_header.split(",")[0], args.out_header.split(",")[1]]
df_q.to_csv(args.out, sep="\t", index=False)


