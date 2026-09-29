#!/usr/bin/env bash

repeats=../../data/2026_09-repeats_grch37/rmsk.txt.gz
genes=../2026_09-living_bed/grch37.genes.sort.bed

zcat $repeats | cut -f 6-8 | sed 's|chr||' > repeat.bed
bedtools coverage -a $genes -b repeat.bed > gene_repeat.coverage.bed

./hist_repeat_cov.py
