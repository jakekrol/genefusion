#!/usr/bin/env bash

gtf=../../data/2026_06-grch37-transcripts/gencode.v19.annotation.gtf.gz

echo "# decompressing gtf"
gunzip -c $gtf > gencode.v19.annotation.gtf

echo "# removing header"
grep -v "^#" gencode.v19.annotation.gtf > gencode.v19.annotation.gtf.nohead

echo "# filtering for exons"
awk -F'\t' 'BEGIN {OFS="\t"} $3 == "exon" {print $1, $4, $5, $7, $9 ,$10}' gencode.v19.annotation.gtf.nohead > gencode.v19.annotation.gtf.nohead.exon

echo "# extracting gene names, trasncript ids, and exon ids"
# produces gencode.v19.annotation.gtf.exons.bed
./extract_gene_name.py \
	--input gencode.v19.annotation.gtf.nohead.exon \
	--output gencode.v19.annotation.gtf.exons.bed \
	--extract_gene_name \
	--extract_transcript_id \
	--extract_exon_id

echo "# sorting exons by chrom, start, end"
sort -t$'\t' -k1,1 -k2,2n -k3,3n gencode.v19.annotation.gtf.exons.bed -o gencode.v19.annotation.gtf.exons.bed.sorted

echo "# compressing exons"
gzip -k gencode.v19.annotation.gtf.exons.bed.sorted

# do the same but for gene-level only

echo "# filtering for genes"
awk -F'\t' 'BEGIN {OFS="\t"} $3 == "gene" {print $1, $4, $5, $7, $9 ,$10}' gencode.v19.annotation.gtf.nohead > gencode.v19.annotation.gtf.nohead.gene

echo "# extracting gene names"
./extract_gene_name.py \
	--input gencode.v19.annotation.gtf.nohead.gene \
	--output gencode.v19.annotation.gtf.nohead.gene.bed \
	--extract_gene_name \
	--extract_transcript_id

sort -t$'\t' -k1,1 -k2,2n -k3,3n gencode.v19.annotation.gtf.nohead.gene.bed -o gencode.v19.annotation.gtf.nohead.gene.bed.sorted

gzip -k gencode.v19.annotation.gtf.nohead.gene.bed.sorted
