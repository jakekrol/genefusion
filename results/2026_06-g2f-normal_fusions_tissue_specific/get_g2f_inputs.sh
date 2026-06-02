#!/usr/bin/env bash

bed='../2026_06-normal_fusions_tissue_specific_set/grch37.genes.added.bed'
fusions='../2026_06-normal_fusions_tissue_specific_set/recurrent_normal_fusion_set.filtered.nohead.c1c2.tsv'
shardfile='../2026_05-g2f-tumor_and_normal-rerun/shardfile.tsv'
cp $bed $fusions $shardfile .