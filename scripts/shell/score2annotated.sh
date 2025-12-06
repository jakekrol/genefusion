#!/usr/bin/env bash
set -euo pipefail

t_0=$(date +%s)
GENOME_LIB_DIR="/data/jake/FusionAnnotator/genome_lib_dir"
echo "# parsing arguments"
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
        -t|--tempdir)
            temp_dir="$2"
            export TMPDIR="$temp_dir"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done
tmp=$(mktemp fusion_annotator_tmp.XXXXXX) || { echo "# failed to create temporary file"; exit 1; }
trap 'rm -f "$tmp"' EXIT


if [ "$header" = true ]; then
    echo "# processing with header"
    # make fusion key column for fusionannnotator (c1--c2)
    echo "# adding a fusion key column"
    awk -v col1="$column1" -v col2="$column2" \
        'BEGIN{FS=OFS="\t"} NR==1{print "FusionKey", $0} NR>1{print $col1"--"$col2, $0}' $input_file > $tmp
    echo "# running FusionAnnotator"
    FusionAnnotator --genome_lib_dir $GENOME_LIB_DIR \
        --fusion_name_col 0 \
        --annotate $tmp \
        > $output_file
else
    echo "# processing without header"
    # make fusion key column for fusionannnotator (c1--c2)
    echo "# adding a fusion key column"
    awk -v col1="$column1" -v col2="$column2" \
        'BEGIN{FS=OFS="\t"} {print $col1"--"$col2, $0}' $input_file > $tmp
    echo "# running FusionAnnotator"
    FusionAnnotator --genome_lib_dir $GENOME_LIB_DIR \
        --fusion_name_col 0 \
        --no_add_header_column \
        --annotate $tmp \
        > $output_file
fi
echo "# annotation complete. output saved to $output_file"

t_1=$(date +%s)
elapsed=$((t_1 - t_0))
echo "# elapsed time: $elapsed seconds"
