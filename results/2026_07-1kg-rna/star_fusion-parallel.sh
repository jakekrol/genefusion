#!/usr/bin/env bash

# default threads for star is 4
PROCESSES=5
export GENOME_LIB_DIR=/data/jake/FusionAnnotator/genome_lib_dir
export OUTDIR=star_fusion-parallel
mkdir -p $OUTDIR

# source ~/.bashrc to init conda
tail -n +2 star_fusion_queries.tsv |
    gargs --log=star_fusion.log \
        -p $PROCESSES \
        "source ~/.bashrc && conda activate star_fusion && STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR --left_fq {1} --right_fq {2} --output_dir $OUTDIR/{0}"

### single example
# GENOME_LIB_DIR=/data/jake/FusionAnnotator/genome_lib_dir
# OUTDIR=star_fusion-SRR19762177
# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
#     --left_fq SRR19762177_1.fastq.gz \
#     --right_fq SRR19762177_2.fastq.gz \
#     --output_dir $OUTDIR