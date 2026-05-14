#!/usr/bin/env bash

# fasta
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# gtf
wget https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gunzip Homo_sapiens.GRCh37.87.gtf.gz

# index
STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir grch37_index \
--genomeFastaFiles Homo_sapiens.GRCh37.dna.primary_assembly.fa \
--sjdbGTFfile Homo_sapiens.GRCh37.87.gtf \
--sjdbOverhang 99 # max.read.len - 1