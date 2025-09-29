#!/usr/bin/env bash

# get queries
bed=../../results/2025_09-gene_bed/grch37.bed
cut -f4 $bed | sort | uniq > queries.txt

# apply gene2tissue.py
./gene2tissue.py \
    -i queries.txt \
    -d ../../data/2025_10-human_protein_atlas_gene_rna_expr/rna_tissue_consensus.tsv \
    -a ../2025_09-haas_interval_frac/merged_alias.yaml \
    -o gene_tissue_top3.tsv \
    -k 3 
