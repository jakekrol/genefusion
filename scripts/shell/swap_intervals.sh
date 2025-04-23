#!/usr/bin/env bash

# input: giggle BED file
# output: BED file with left and right intervals swapped (c1-4 -> c5-8 & c5-8 -> c1-4)

in=$1
out=$2
echo "Swapping intervals in $in to $out"
awk 'BEGIN {OFS="\t"} {print $5, $6, $7, $8, $1, $2, $3, $4, $9, $10}' $in > $out
