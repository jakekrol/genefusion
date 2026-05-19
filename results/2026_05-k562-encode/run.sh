#!/usr/bin/env bash

# build grch37 star index for alignment
./star_index.sh
# align with star and get chimeric reads
./star_align.sh
# convert chimeric reads to bedpe format
./star2bedpe.py
# index
mkdir -p k562-rna-bed
cp k562-rna-fusion.bedpe.gz k562-rna-bed/k562.bed.gz
giggle index -s -i "k562-rna-bed/*.bed.gz" -o k562-rna-index

