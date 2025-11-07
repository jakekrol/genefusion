#!/usr/bin/env bash

# input: swapped giggle BED file
# output: intersected BED file with gene_file.txt.latest

export bedtools='/data/jake/bedtools.static.binary'
giggleswap=$1
genebed=$2
out=$3

echo "Intersecting $giggleswap with $genebed to $out"
"$bedtools" intersect -a $genebed \
    -b $giggleswap -wb -wa > "$out"

