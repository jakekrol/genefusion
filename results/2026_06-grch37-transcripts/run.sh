#!/usr/bin/env bash

gtf=../../data/2026_06-grch37-transcripts/gencode.v19.annotation.gtf.gz

gunzip -c $gtf > gencode.v19.annotation.gtf

grep -v "^#" gencode.v19.annotation.gtf > gencode.v19.annotation.gtf.nohead

awk -F'\t' 'BEGIN {OFS="\t"} $3 == "exon" {print $1, $4, $5, $7, $9 ,$10}' gencode.v19.annotation.gtf.nohead > gencode.v19.annotation.gtf.nohead.exon

# produces gencode.v19.annotation.gtf.exons.bed
./extract_gene_name.py

gzip -k gencode.v19.annotation.gtf.exons.bed




