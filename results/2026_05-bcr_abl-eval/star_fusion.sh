#!/usr/bin/env bash

GENOME_LIB_DIR="/data/jake/FusionAnnotator/genome_lib_dir"

fractions=(0.2 0.4 0.6 0.8)
for fraction in "${fractions[@]}"; do
    STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
                 --left_fq "k562-${fraction}-1.fastq" \
                 --right_fq "k562-${fraction}-2.fastq" \
                 --output_dir "star_fusion_outdir_${fraction}"
done