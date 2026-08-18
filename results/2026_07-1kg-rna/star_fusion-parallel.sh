#!/usr/bin/env bash

# ! conda activate star_fusion

# default threads for star is 4
CPUS=10
export GENOME_LIB_DIR=/data/jake/FusionAnnotator/genome_lib_dir
export OUTDIR=star_fusion-parallel

tail -n +2 star_fusion_queries.tsv |
    gargs --log=star_fusion.log \
        -p $CPUS \
        "STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR --left_fq {1} --right_fq {2} --output_dir $OUTDIR/{0}"

GENOME_LIB_DIR=/data/jake/FusionAnnotator/genome_lib_dir
OUTDIR=star_fusion-SRR19762177
STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
    --left_fq SRR19762177_1.fastq.gz \
    --right_fq SRR19762177_2.fastq.gz \
    --output_dir $OUTDIR