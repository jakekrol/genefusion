#!/usr/bin/env bash
CPUS=16

STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
			--left_fq SRR19762165_1.fastq.gz \
			--right_fq SRR19762165_2.fastq.gz \
			--output_dir SRR19762165-star_fusion-out \
            --CPU $CPUS

samtools sort -@ $CPUS -o SRR19762165.star_fusion.sort.bam SRR19762165-star_fusion-out/Aligned.out.bam
samtools index SRR19762165.star_fusion.sort.bam