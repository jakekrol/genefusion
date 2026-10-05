#!/usr/bin/env bash

GENOME_LIB_DIR=$GENOME_LIB_DIR
mapfile -t files < <(ls top_fusions/*)
rm annotate.input || echo "# no annotate.input file found. creating new."
for x in "${files[@]}"; do
    echo $x >> annotate.input
done

for f in "${files[@]}"; do
    out="${f%.tsv}.annotated.tsv"
    tail -n +2 $f | \
        awk -v OFS="\t" '{print $1"--"$2, $3}' | \
            FusionAnnotator --genome_lib_dir $GENOME_LIB_DIR --no_add_header_column | \
                grep -vE "SELFIE|NEIGHBOR|BLAST|PARALOGS" > \
                $out
done
