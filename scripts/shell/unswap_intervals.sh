#!/usr/bin/env bash

# input: intersect BED file with 
# - c1-5 is from genefile
#         - c6-9 is rh
#         -  c10-13 is lh
# - c14 is ?
# - c15 is sample
# output: intersect BED file with lh at c6-9 and rh at c10-13

in=$1
out=$2
echo "Unswapping intervals from $in to $out"
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $10, $11, $12, $13, $6, $7, $8, $9, $14, $15}' $in > $out
