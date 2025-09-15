#!/usr/bin/env bash

# setup haas for join with 1kg results
./setup_haas.py \
    --haas ../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt \
    --bedfile ../results/2025_09-gene_bed/grch37.bed \
    --alias_map ../2025_09-haas_interval_frac/merged_alias.yaml \
    --output 'haas_left_right_mapped.tsv'

# join haas with 1kg results
onekg="../results/2025_09-onekg_giggle2fusion/g2f/pop_normal_fusions.tsv"

join.py \
    -x haas_left_right_mapped.tsv \
    -y $onekg \
    -t left \
    -o haas_onekg.tsv \
    -k left,right



