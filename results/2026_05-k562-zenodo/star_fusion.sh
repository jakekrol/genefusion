#!/usr/bin/env bash

GENOME_LIB_DIR="/data/jake/FusionAnnotator/genome_lib_dir"
# STAR-Fusion --genome_lib_dir "${GENOME_LIB_DIR}" \
#             -J k562-rna-fusion-Chimeric.out.junction \
#             --output_dir star_fusion_outdir

 STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
             --left_fq k562_1.fastq \
             --right_fq k562_2.fastq \
             --output_dir star_fusion_outdir2
