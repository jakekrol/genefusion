#!/usr/bin/env bash

dir_chunks_scored='chunks_scored'
mapfile -t tissues < <(ls "$dir_chunks_scored")
for tissue in "${tissues[@]}"; do
	echo "# tissue: $tissue"
	outfile="$dir_chunks_scored/${tissue}_scored.tsv"
	echo "# outfile: $outfile"
	printf "gene_left\tgene_right\tscore_uniform\n" > "$outfile"
	dir_tissue="$dir_chunks_scored/$tissue"
	mapfile -t chunks < <(ls "$dir_tissue")
	for chunk in "${chunks[@]}"; do
		tail -n +2 "$dir_tissue/$chunk" >> "$outfile"
	done
	sort_tbl.sh -i "$outfile" -o "$outfile.sorted" -c score_uniform -s DESC
done