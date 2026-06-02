#!/usr/bin/env bash

join_script='join.py'
fout='all-evidence.tsv'
mapfile -t files < <(ls g2f_agg/*.tsv)


f1="${files[0]}"
cp "$f1" "$fout"
files=("${files[@]:1}")
for f in "${files[@]}"; do
    $join_script \
        --type outer \
        --keys gene_left,gene_right \
        -x "$fout" \
        -y "$f" \
        -o "$fout" \
        --verbose
done

# keep only relevant fusions

join.py \
    --type left \
    --keys gene_left,gene_right \
    -x query_fusions_sorted.tsv \
    -y "$fout" \
    -o query-evidence.tsv
    --verbose

# fillna

fillna.sh \
    -i query-evidence.tsv \
    -o query-evidence-filled.tsv \
    -f 0