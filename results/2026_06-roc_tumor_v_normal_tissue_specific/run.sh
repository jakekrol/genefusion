#!/usr/bin/env bash

script=roc.py

normal_agg_subpop_weighted=../2026_06-g2f-normal_fusions_tissue_specific/normal_score.agg.subpop.weighted.tsv
normal_agg_subpop_uniform=../2026_06-g2f-normal_fusions_tissue_specific/normal_score.agg.subpop.uniform.tsv
normal_in_tissue=../2026_06-g2f-normal_fusions_tissue_specific/normal_score.in_tissue.tsv

tumor_agg_subpop_weighted=../2026_05-g2f-tumor_and_normal-rerun/tumor_score.agg.subpop.weighted.tsv
tumor_agg_subpop_uniform=../2026_05-g2f-tumor_and_normal-rerun/tumor_score.agg.subpop.uniform.tsv
tumor_in_tissue=../2026_05-g2f-tumor_and_normal-rerun/tumor_score.in_tissue.tsv

tumor_files=("${tumor_agg_subpop_weighted}" "${tumor_agg_subpop_uniform}" "${tumor_in_tissue}")
for f in "${tumor_files[@]}"; do
	fbase="$(basename $f)"
	printf "fusion\tscore\tlabel\n" > $fbase
	tail -n +2 $f | sed 's|$|\t1|g' >> $fbase
done

normal_files=("${normal_agg_subpop_weighted}" "${normal_agg_subpop_uniform}" "${normal_in_tissue}")
for f in "${normal_files[@]}"; do
	fbase="$(basename $f)"
	printf "fusion\tscore\tlabel\n" > $fbase
	tail -n +2 $f | sed 's|$|\t0|g' >> $fbase
done

cat tumor_score.agg.subpop.weighted.tsv \
	<(tail -n +2 normal_score.agg.subpop.weighted.tsv) > score.agg.subpop.weighted.tsv
cat tumor_score.agg.subpop.uniform.tsv \
	<(tail -n +2 normal_score.agg.subpop.uniform.tsv) > score.agg.subpop.uniform.tsv
cat tumor_score.in_tissue.tsv \
	<(tail -n +2 normal_score.in_tissue.tsv) > score.in_tissue.tsv

roc.py \
	--scores "score.agg.subpop.weighted.tsv,score.agg.subpop.uniform.tsv,score.in_tissue.tsv" \
	--names "agg_subpop_weighted,agg_subpop_uniform,in_tissue" \
	--header \
	--score_col 1 \
	--label_col 2 \
	--output roc.png \
	--reference

### scatter

# join evidence and score
tumor_evidence=../2026_05-g2f-tumor_and_normal-rerun/query-evidence-filled.tsv
normal_evidence=../2026_06-g2f-normal_fusions_tissue_specific/query-evidence-filled.tsv

cat $tumor_evidence <(tail -n +2 $normal_evidence ) > all_evidence.tsv
script="import pandas as pd
df=pd.read_csv('all_evidence.tsv', sep='\t')
df['fusion'] = df['gene_left'] + '--' + df['gene_right']
df.to_csv('all_evidence.tsv', sep='\t', index=False)
"
python -c "$script"
join.py \
	-x score.in_tissue.tsv \
	-y all_evidence.tsv \
	--keys fusion \
	--type left \
	--output all_evidence_with_score.tsv

# join tissues


