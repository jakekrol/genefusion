#!/usr/bin/env bash

# wgs
# >100GB
# paired-end, 151nt
wget 'https://www.encodeproject.org/files/ENCFF287PIC/@@download/ENCFF287PIC.bam'
samtools index ENCFF287PIC.bam

# rna-seq
# ~4 GB each
# 2 runs with same library
# single-end, 100nt
wget 'https://www.encodeproject.org/files/ENCFF764ZJT/@@download/ENCFF764ZJT.fastq.gz'
wget 'https://www.encodeproject.org/files/ENCFF133KON/@@download/ENCFF133KON.fastq.gz'