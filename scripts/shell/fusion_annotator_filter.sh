#!/usr/bin/env bash
set -euo pipefail

GENOME_LIB_DIR="/data/jake/FusionAnnotator/genome_lib_dir"
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input_file="$2"
            shift 2
            ;;
        -o|--output)
            output_file="$2"
            shift 2
            ;;
        --header)
            header=true
            shift 1
            ;;
        -c1|--column1)
            column1="$2"
            shift 2
            ;;
        -c2|--column2)
            column2="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done


if [ "$header" = true ]; then
    # make fusion key column for fusionannnotator (c1--c2)
    awk -v col1="$column1" -v col2="$column2" \
        'BEGIN{FS=OFS="\t"} NR==1{print "FusionKey", $0} NR>1{print $col1"--"$col2, $0}' $input_file > temp_input.txt
    FusionAnnotator --genome_lib_dir $GENOME_LIB_DIR \
        --fusion_name_col 0 \
        --annotate temp_input.txt | grep -v PARA | grep -v NEIGHBORS_OVERLAP | grep -v BLASTPAIR | grep -v Normal \
        > $output_file
else
    awk -v col1="$column1" -v col2="$column2" \
        'BEGIN{FS=OFS="\t"} {print $col1"--"$col2, $0}' $input_file > temp_input.txt
    FusionAnnotator --genome_lib_dir $GENOME_LIB_DIR \
        --fusion_name_col 0 \
        --no_add_header_column \
        --annotate temp_input.txt | grep -v PARA | grep -v NEIGHBORS_OVERLAP | grep -v BLASTPAIR | grep -v Normal \
        > $output_file
fi
rm temp_input.txt
echo "Annotation complete. Output saved to $output_file"
