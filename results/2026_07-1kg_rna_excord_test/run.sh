#!/usr/bin/env bash

sort_bed='/data/jake/genefusion/scripts/shell/sort_bed'
CPUS=16

./download.py

# conda activate star_fusion
./star_fusion.sh

samtools sort -@ $CPUS -o SRR19762165.star_fusion.sort.bam SRR19762165-star_fusion-out/Aligned.out.bam
samtools index SRR19762165.star_fusion.sort.bam

./excord.sh

# sort
mkdir tmp1 tmp2
mv SRR19762165.star_fusion.sort.excord.bed.gz tmp1
$sort_bed tmp1 tmp2 $CPUS

# index
mkdir -p bed_sort
mv tmp2/* bed_sort
giggle index -s -i "bed_sort/*.bed.gz" -o giggle_index

# search
gene_left=TIMM23
gene_left_chr=10
gene_left_start=51592080
gene_left_end=51623365
region=$gene_left_chr:$gene_left_start-$gene_left_end
giggle search -i giggle_index -r $region -v > $gene_left.giggle

# intersect
bedtools intersect -a $GRCH37_GENES -b $gene_left.giggle -wa -wb > $gene_left.intersect

# count genes
cut -f 4 $gene_left.intersect | sort | uniq -c | sort -nr > $gene_left.counts




