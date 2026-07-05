#!/usr/bin/env bash

CPU=25
GENOME_LIB_DIR=/data/jake/FusionAnnotator/genome_lib_dir

# # kpl4
# kpl4_1='../../data/2026_06-cancer_cell_line-fastq/kpl4_1.fastq.gz'
# kpl4_2='../../data/2026_06-cancer_cell_line-fastq/kpl4_2.fastq.gz'

# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq $kpl4_1 \
# 			 --right_fq $kpl4_2 \
# 			 --output_dir star_fusion_outdir-kpl4 \
# 			 --CPU ${CPU}

# # mcf7
# mcf7_1='../../data/2026_06-cancer_cell_line-fastq/mcf7_1.fastq.gz'
# mcf7_2='../../data/2026_06-cancer_cell_line-fastq/mcf7_2.fastq.gz'

# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq $mcf7_1 \
# 			 --right_fq $mcf7_2 \
# 			 --output_dir star_fusion_outdir-mcf7 \
# 			 --CPU ${CPU}

# # vcap85
# vcap85_1='../../data/2026_06-cancer_cell_line-fastq/vcap85_1.fastq.gz'
# vcap85_2='../../data/2026_06-cancer_cell_line-fastq/vcap85_2.fastq.gz'

# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq $vcap85_1 \
# 			 --right_fq $vcap85_2 \
# 			 --output_dir star_fusion_outdir-vcap85 \
# 			 --CPU ${CPU}

# # lc2ad
# lc2ad_1='../../data/2026_06-cancer_cell_line-fastq/lc2ad_1.fastq.gz'
# lc2ad_2='../../data/2026_06-cancer_cell_line-fastq/lc2ad_2.fastq.gz'

# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq $lc2ad_1 \
# 			 --right_fq $lc2ad_2 \
# 			 --output_dir star_fusion_outdir-lc2ad \
# 			 --CPU ${CPU}

# # h2228
# h2228_1='../../data/2026_06-cancer_cell_line-fastq/h2228_1.fastq.gz'
# h2228_2='../../data/2026_06-cancer_cell_line-fastq/h2228_2.fastq.gz'

# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq $h2228_1 \
# 			 --right_fq $h2228_2 \
# 			 --output_dir star_fusion_outdir-h2228 \
# 			 --CPU ${CPU}

# # skbr3
# skbr3_1='../../data/2026_06-cancer_cell_line-fastq/skbr3_1.fastq.gz'
# skbr3_2='../../data/2026_06-cancer_cell_line-fastq/skbr3_2.fastq.gz'

# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq $skbr3_1 \
# 			 --right_fq $skbr3_2 \
# 			 --output_dir star_fusion_outdir-skbr3 \
# 			 --CPU ${CPU}

# bt474
bt474_1='../../data/2026_06-cancer_cell_line-fastq/BT474_1.fastq.gz'
bt474_2='../../data/2026_06-cancer_cell_line-fastq/BT474_2.fastq.gz'

STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
			 --left_fq $bt474_1 \
			 --right_fq $bt474_2 \
			 --output_dir star_fusion_outdir-bt474 \
			 --CPU ${CPU}