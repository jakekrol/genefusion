#!/usr/bin/env bash

# setup haas for join with 1kg results
./setup_haas.py \
    --haas ../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt \
    --bedfile ../results/2025_09-gene_bed/grch37.bed \
    --alias_map ../2025_09-haas_interval_frac/merged_alias.yaml \
    --output 'haas_left_right_mapped.tsv'

# rm unresolved gene names
grep -v "^-1" haas_left_right_mapped.tsv > haas_left_right_mapped_cln.tsv

# join haas with 1kg results
onekg="../2025_09-onekg_giggle2fusion/g2f/pop_normal_fusions.tsv"

join.py \
    -x haas_left_right_mapped_cln.tsv \
    -y $onekg \
    -t left \
    -o haas_onekg.tsv \
    -k left,right

drop_duplicates.py -i haas_onekg.tsv -o haas_onekg_no_dups.tsv -c 5

# histogram of 1kg supporting reads for haas fusions

tail -n +2 haas_onekg_no_dups.tsv \
    | cut -f11 \
    | hist.py -o haas_onekg_readcount_histogram.png \
    --ylog -x "1000 Genomes depth" -y Frequency

# inspect fusions with >10000 supporting reads in 1kg
tail -n +2 haas_onekg_no_dups.tsv \
    | awk '$11 > 10000' | \
    # nr numeric reverse
    sort -k11,11nr | \
    # reorder columns to put depth third
    awk 'BEGIN{OFS="\t"} {print $1,$2,$11,$3,$4,$5,$6,$7,$8,$9,$10}' > haas_onekg_depth_gt10000.tsv

# do the same with >1000 supporting reads in 1kg
tail -n +2 haas_onekg_no_dups.tsv \
    | awk '$11 > 1000' | \
    sort -k11,11nr | \
    awk 'BEGIN{OFS="\t"} {print $1,$2,$11,$3,$4,$5,$6,$7,$8,$9,$10}' > haas_onekg_depth_gt1000.tsv  
