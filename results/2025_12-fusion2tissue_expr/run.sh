#!/usr/bin/env bash

set -euo pipefail
t_0=$(date +%s)
K=3
fusion_key=../2025_11-fusion_lexicographic_key/gene_pairs_sorted_by_genomic_position.tsv

echo "# using K=$K for top K tissues"
echo "# using fusion key: $fusion_key"

# echo "# getting unique genes from fusion key"
# tmp=$(mktemp "fusion2tissue_XXXX.tsv")
# trap 'rm -f "$tmp"' EXIT
# echo "# tmp file: $tmp"
# tail -n +2 "$fusion_key" | cut -f1 > "$tmp"
# tail -n +2 "$fusion_key" | cut -f2 >> "$tmp"
# sort -u "$tmp" > "./genes_unique.tsv"

echo "# mapping genes to top tissue expression"
../../scripts/python/gene2tissue.py \
    --input ./genes_unique.tsv \
    --data ../../data/2025_10-human_protein_atlas_gene_rna_expr/rna_tissue_consensus.tsv \
    --topk $K \
    --alias ../2025_09-haas_interval_frac/merged_alias.yaml \
    --out ./gene2tissue_top${K}_expr.tsv \
    --out_header gene,tissues_top${K}_expr

t_1=$(date +%s)
echo "# done. total time: $((t_1 - t_0)) seconds"





