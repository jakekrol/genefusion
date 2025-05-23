#!/usr/bin/env bash

# input: unswapped intersected GIGGLE BED file
# output: 2-column file: c1 is fusion gene, c2 is PE evidence count

echo "Deprecated: use count_fusions.py instead"
in=$1
out=$2
cut -f4 $in | sort | uniq -c  |  awk 'BEGIN {OFS="\t"} {print $2, $1}'  > $out
