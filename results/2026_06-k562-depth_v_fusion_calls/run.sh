#!/usr/bin/env bash

SEED=1
GENOME_LIB_DIR="/data/jake/FusionAnnotator/genome_lib_dir"
CPU=25
OUTFILE="fusion_counts.tsv"


# echo "# gathering and decompressing fastq files..."
# fastq1='../../data/2026_05-k562-ccle/SRR521460_1.fastq.gz'
# fastq2='../../data/2026_05-k562-ccle/SRR521460_2.fastq.gz'
# gunzip -k -c "${fastq1}" > k562_1.fastq
# gunzip -k -c "${fastq2}" > k562_2.fastq

# fractions=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
# echo "downsampling fastq files to fractions: ${fractions[@]}"
# for fraction in "${fractions[@]}"; do
# 	echo "# downsampling fastq files to fraction ${fraction}..."
# 	seqtk sample -s${SEED} k562_1.fastq ${fraction} > k562_1-${fraction}.fastq
# 	seqtk sample -s${SEED} k562_2.fastq ${fraction} > k562_2-${fraction}.fastq
# done

# echo "# running STAR-Fusion on downsampled fastq files..."
# fractions=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
# for fraction in "${fractions[@]}"; do
# 	echo "# running STAR-Fusion on downsampled fastq files at fraction ${fraction}..."
# 	STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq k562_1-${fraction}.fastq \
# 			 --right_fq k562_2-${fraction}.fastq \
# 			 --output_dir star_fusion_outdir-${fraction} \
# 			 --CPU ${CPU}
# done
# echo "# run on full fastq files..."
# STAR-Fusion --genome_lib_dir $GENOME_LIB_DIR \
# 			 --left_fq k562_1.fastq \
# 			 --right_fq k562_2.fastq \
# 			 --output_dir star_fusion_outdir-1.0 \
# 			 --CPU ${CPU}


echo "# counting number of fusion called per fraction..."
fractions=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
rm ${OUTFILE} || echo "# ${OUTFILE} does not exist, creating new file..."
for fraction in "${fractions[@]}"; do
	count=$(grep -v \
		"^#" star_fusion_outdir-${fraction}/star-fusion.fusion_predictions.abridged.tsv \
		| wc -l)
	printf "${fraction}\t${count}\n" >> ${OUTFILE}
done

echo "# plotting depth versus fusion counts..."
./plot.py
