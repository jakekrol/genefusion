#!/usr/bin/env bash

../../scripts/python/build_gene_bed.py \
    --hgnc ../../data/2025_09-gene_bed/genenames.tsv \
    --gtf ../../data/2025_09-gene_bed/gencode.v19.annotation.gtf.gz \
    --output grch37.chr.bed

# rm mitochondrial genes
grep -v '^chrM' grch37.chr.bed > tmp && mv tmp grch37.chr.bed

# make a copy with "chr" removed from chromosome names
sed 's|^chr||' grch37.chr.bed > grch37.bed

# rm '.' column
cut -f 5 --complement grch37.bed > tmp && mv tmp grch37.bed

# replace '+/-' with pos/neg
sed -i 's|+|pos|g; s|-|neg|g' grch37.bed