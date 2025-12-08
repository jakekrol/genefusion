#!/usr/bin/env bash

set -euo pipefail
t_0=$(date +%s)
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
        --gene2tissue_topk)
            gene2tissue_topk="$2"   
            shift 2
            ;;
        --tempdir)
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

echo "# using tmpdir: $TMPDIR"
tmp_g2t=$(mktemp "$TMPDIR/score2tissue_expr_tmp_topk_file.XXXXXX") || { echo "# failed to create temporary file"; exit 1; }
tmp_left=$(mktemp "$TMPDIR/score2tissue_expr_tmp_left_file.XXXXXX") || { echo "# failed to create temporary file"; exit 1; }
trap 'rm -f "$tmp_g2t"' EXIT
trap 'rm -f "$tmp_left"' EXIT
echo "# getting tissue for left gene from $gene2tissue_topk"
echo "# using temporary file: $tmp_g2t"
# format the header
sed '1s/.*/left\tleft_top_expr/' "$gene2tissue_topk" > "$tmp_g2t"
# join left gene with gene2tissue_topk
join_ddb.py \
    --left "$input_file" \
    --right "$tmp_g2t" \
    --type left \
    --key left \
    --output "$tmp_left"
echo "# getting tissue for right gene from $gene2tissue_topk"
# format the header
sed '1s/.*/right\tright_top_expr/' "$gene2tissue_topk" > "$tmp_g2t"
# join right gene with gene2tissue_topk
join_ddb.py \
    --left "$tmp_left" \
    --right "$tmp_g2t" \
    --type left \
    --key right \
    --output "$output_file"

t_1=$(date +%s)
echo "# done. total time: $((t_1 - t_0)) seconds"
