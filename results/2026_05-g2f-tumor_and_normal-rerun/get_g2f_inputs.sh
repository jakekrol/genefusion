#!/usr/bin/env bash

# shardfile
shardfile='../2026_04-s2f_babiceanu_recurrent_normal_tissue_agnostic/shardfile.tsv'
cut -f 1,4 $shardfile > z && mv z shardfile.tsv

# normal fusion set
normal_fusions='../2026_04-s2f_babiceanu_recurrent_normal_tissue_agnostic/recurrent_normal_fusions.tsv'
cp $normal_fusions .
normal_fusions='recurrent_normal_fusions.tsv'

# tumor fusion set
tumor_fusions='../2026_04-s2f_pcawg_recurrent/recurrent_tumor_fusions_set.tsv'
cp $tumor_fusions .
tumor_fusions='recurrent_tumor_fusions_set.tsv'

cat $normal_fusions $tumor_fusions > query_fusions.tsv

# bedfile
bedfile='../2026_04-s2f_babiceanu_recurrent_normal_tissue_agnostic/grch37.genes.bed.added'
cp $bedfile .