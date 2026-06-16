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
