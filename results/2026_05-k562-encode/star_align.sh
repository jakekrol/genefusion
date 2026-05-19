#!/usr/bin/env bash

# fastq1='../../data/2026_05-k562/ENCFF133KON.fastq.gz'
# fastq2='../../data/2026_05-k562/ENCFF764ZJT.fastq.gz'
# gunzip -k -c "${fastq1}" > ENCFF133KON.fastq
# gunzip -k -c "${fastq2}" > ENCFF764ZJT.fastq
fastq1='ENCFF133KON.fastq'
fastq2='ENCFF764ZJT.fastq'
echo "# aligning ${fastq1} and ${fastq2} with STAR..."
genome_index='./grch37_index'
# info on how to run star
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
# info on how to run star to include necessary files for fusion detection
# https://github.com/STAR-Fusion/STAR-Fusion/wiki#alternatively-kickstart-mode-running-star-yourself-and-then-running-star-fusion-using-the-existing-outputs
# when mapping mulitple single-end files jointly use comma-separator
# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf (pg8)
STAR --genomeDir "${genome_index}" \
    --runThreadN 10 \
    --readFilesIn "${fastq1},${fastq2}" \
    --outFileNamePrefix k562-rna-fusion- \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outReadsUnmapped None \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 8 \
    --chimOutJunctionFormat 1 \
    --alignSJDBoverhangMin 10 \
    --alignMatesGapMax 100000 \
    --alignIntronMax 100000 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --outSAMattrRGline ID:GRPundef \
    --chimMultimapScoreRange 3 \
    --chimScoreJunctionNonGTAG -4 \
    --chimMultimapNmax 20 \
    --chimNonchimScoreDropMin 10 \
    --peOverlapNbasesMin 12 \
    --peOverlapMMp 0.1 \
    --alignInsertionFlush Right \
    --alignSplicedMateMapLminOverLmate 0 \
    --alignSplicedMateMapLmin 30

STAR --genomeDir "${genome_index}" \
    --runThreadN 10 \
    --readFilesIn "${fastq1},${fastq2}" \
    --outFileNamePrefix k562-rna- \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterScoreMin 0 \
    --outFilterMatchNmin 0 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 1.0 \
    --outFilterMultimapNmax 999 \
    --outFilterMultimapScoreRange 0