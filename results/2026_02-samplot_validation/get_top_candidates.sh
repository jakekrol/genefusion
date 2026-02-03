#!/usr/bin/env bash

t_0=$(date +%s)

dir_score_base="../2025_11-score"
# tissues_w_rna=(blood esophagus kidney liver ovary) haven't generated final scores for blood and esophagus yet
tissues_w_rna=(kidney liver ovary)
final_score_filename=scored_dna1.0_t0.5_r0.5_u50_expr_anno_sort_no_dups_anno_filt_burden_filt.tsv
k=100

for tissue in "${tissues_w_rna[@]}"; do
    dir_score="${dir_score_base}/${tissue}"
    file_score="${dir_score}/${final_score_filename}"
    outfile="./${tissue}_top_${k}_candidates.tsv"
    echo "# getting top k=${k} candidates for tissue: ${tissue}"
    echo "# score file: ${file_score}"
    echo "# outfile: ${outfile}"
    if [[ ! -f "${file_score}" ]]; then
        echo "Error: score file not found: ${file_score}"
        exit 1
    fi
    grep INTRACHROMOSOMAL "${file_score}" | head -n $k > "${outfile}"
done
t_1=$(date +%s)
echo "# done. total time: $((t_1 - t_0)) sec"