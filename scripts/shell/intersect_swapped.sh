#!/usr/bin/env bash

# input: swapped giggle BED file
# output: intersected BED file with gene_file.txt.latest

export bedtools='/data/jake/bedtools.static.binary'
giggleswap=$1
out=$2

genefile='/data/jake/genefusion/data/gene_file.txt.latest'
echo "Intersecting $giggleswap with $genefile to $out"
"$bedtools" intersect -a $genefile \
    -b $giggleswap -wb -wa > "$out"

