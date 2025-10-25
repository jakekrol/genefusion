#!/usr/bin/env bash

# download gencode v19
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# unzip
gunzip gencode.v19.annotation.gtf.gz


# keep only chr 1-22, X, Y
grep -E '^#|^chr([1-9][0-9]?|X|Y)\s' gencode.v19.annotation.gtf \
    > gencode.v19.annotation.filtered.gtf

# keep genes only
awk '$3 == "gene"' gencode.v19.annotation.filtered.gtf \
    > z && mv z gencode.v19.annotation.filtered.gtf
# rm 
# pseudogenes
# lincRNA
# processed_transcript
# snRNA
# miRNA
# misc_RNA
# snoRNA
# rRNA
# we keep antisense, sense_intronic, sense_overlapping, and protein_coding
grep -v \
    -E "pseudogene|lincRNA|processed_transcript|snRNA|miRNA|misc_RNA|snoRNA|rRNA" \
    gencode.v19.annotation.filtered.gtf > z && mv z gencode.v19.annotation.filtered.gtf

# convert to bed
# make a padded version for promoters, too
script="
import pandas as pd
import sys,os
import swifter
f='gencode.v19.annotation.filtered.gtf'
df = pd.read_csv(f, sep='\t', header=None, usecols=[0,3,4,6,8])
df.columns = ['chrom', 'start', 'end', 'strand', 'attributes']
df['attributes'] = df['attributes'].str.extract('gene_name \"(.*?)\";')
df = df[['chrom', 'start', 'end', 'attributes', 'strand']]
df['strand'] = df['strand'].replace({'+': 'pos', '-': 'neg'})
# rm chr prefix
df['chrom'] = df['chrom'].str.replace('chr', '')
df.to_csv('grch37.genes.bed', sep='\t', header=False, index=False)
# pad for promoters
pad = 1000
def pad_promoter(row):
    if row['strand'] == 'pos':
        # max to avoid negative start
        new_start = max(0, row['start'] - pad)
        row['start'] = new_start
    else:
        new_end = row['end'] + pad
        row['end'] = new_end
    return row
df = df.swifter.apply(pad_promoter, axis=1)
# write padded bed
df.to_csv('grch37.genes.promoter_pad.bed', sep='\t', header=False, index=False)
"
python3 -c "$script"
