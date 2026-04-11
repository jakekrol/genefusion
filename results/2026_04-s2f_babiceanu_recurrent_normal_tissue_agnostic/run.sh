#!/usr/bin/env bash

normal_fusions='../../data/2026_03-babiceanu_normal_fusions/recurrent_normal_tissue_agnostic.clean.csv'

# make a tsv of fusions with gene metadata
cat $normal_fusions | tr ',' '\t' > normal_fusions.tsv

# split fusion column and rm others
tail -n +2 normal_fusions.tsv | cut -f 1 | tr '-' '\t' > t && mv t normal_fusions.tsv

# assign left, right and gene metadata

