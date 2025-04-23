#!/usr/bin/env bash

# note: make sure to rm the dir prefix and other file extensions from sample column before running this script
# sed s"|blood_sort/||" file > filecln
# sed 's/\.excord\.bed\.gz//' filecln > filecln2

# input: unswapped intersected GIGGLE BED file
# output: the same files split by sample column

in=$1
outdir=$2
samplecol=15

in_base=$(basename $in)

echo "splitting $in to sample wise files to $outdir"
# the output file name is sample.input_file_name
awk -F'\t' -v outdir="$outdir" -v suffix="$in_base" -v j=$samplecol '{print > outdir"/"$j"."suffix}' "$in"
