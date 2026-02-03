#!/usr/bin/env bash

./get_top_candidates.sh

for tissue in kidney liver ovary; do
    input_file="./${tissue}_top_100_candidates.tsv"
    output_file="./${tissue}_top_candidates-stix_coords.tsv"
    echo "# processing tissue: ${tissue}"
    echo "# input file: ${input_file}"
    echo "# output file: ${output_file}"
    ./fusion2stix_query.sh -i <(cut -f 1,2 "${input_file}") -o "${output_file}"
done