#!/usr/bin/env bash

tmp='/data/jake/tmp'
bedfile='../2025_04-gene_bedfile_cln/grch37.genes.bed'
# find fusion with highest read depth in 1000G high coverage
fusion_evidence_1000g='../2026_05-g2f-onekg_high_cov/g2f_agg/high_coverage_1000g_dna-fusion_evidence.tsv'

# filter out selfies
awk -F'\t' -v OFS='\t' 'NR>1 && $1 != $2 { print }' $fusion_evidence_1000g | \
    LC_ALL=C sort -k 3,3nr --buffer-size 3G -T $tmp | \
    head -n 50 > highest_read_depth_fusion.tsv

cut -f 1,2 highest_read_depth_fusion.tsv > highest_read_depth_fusion_genepairs.tsv

# find any without overlapping annotations
./find_non_overlapping_annotations.py

# from here, i found NLK and PIGS
gene_left=NLK
gene_right=PIGS
region_nlk=17:26368763-26523407
region_pigs=17:26880401-26898890

#get the evidence for these two genes
nlk_intersect=../2026_05-g2f-onekg_high_cov/g2f_out/high_coverage_1000g_dna/NLK.giggle.clean.swap.intersect.bed.gz
printf "gene_chr_right\tgene_start_right\tgene_end_right\tgene_name_right\tgene_strand_right\tchr_right\tstart_right\tend_right\tstrand_right\tchr_left\tstart_left\tend_left\tstrand_left\tevidence_type\tsample\n" > nlk_pigs_intersect.bed
zcat $nlk_intersect | grep -w $gene_right >> nlk_pigs_intersect.bed




