#!/usr/bin/env bash

# conda activate polymerization

set -euo pipefail

CPUS=30
JOIN_SCRIPT=/data/jake/rl-tools/wrangle/join_ddb.py
FINAL_TABLE=all_categories-fusion_table.tsv

recurrent_tumor_fusions='../../data/2026_04-pcawg_recurrent_fusions/recurrent_tumor_fusions.csv'
# csv to tsv
python -c "import pandas as pd; df = pd.read_csv('$recurrent_tumor_fusions'); df.to_csv('recurrent_tumor_fusions.tsv', sep='\t', index=False)"
# make fusion set
tail -n +2 recurrent_tumor_fusions.tsv | cut -f 1,2 | sort -u > recurrent_tumor_fusions_set.tsv
# get stix shards
cp ../../config/shardfile.tsv .
# get bedfile
cp ../2025_04-gene_bedfile_cln/grch37.genes.bed .

./run_s2f.py \
	--shardfile shardfile.tsv \
	--fusionset recurrent_tumor_fusions_set.tsv \
	--bedfile grch37.genes.bed \
	--cpus $CPUS \
	--outdir stix_output

# aggregate by index category
./agg_within_index_category.py 2>&1 | tee agg_within_index_category.log
# join all results into one file
shopt -s nullglob
fusiontables=(stix_output_agg/*.tsv)

if ((${#fusiontables[@]} == 0)); then
	echo "No TSV files found in stix_output_agg/" >&2
 	exit 1
fi

# select first table to join on, then iteratively join the rest
tbl1=${fusiontables[0]}
cp -- "$tbl1" "$FINAL_TABLE"
for tbl in "${fusiontables[@]:1}"; do
	# iteration 1 uses first table as left
	# iteration >= 2 uses final output table path as left
	# right table is always the table selected by iteration
	echo "# joining $FINAL_TABLE and $tbl"
 	"$JOIN_SCRIPT" \
		--left "$FINAL_TABLE" \
		--right "$tbl" \
		--keys gene_left,gene_right \
		--type outer \
		--output "$FINAL_TABLE"
done

# get score column maps
cp ../../config/tumor_colmap.yaml ../../config/normal_colmap.yaml .
./run_score.py 2>&1 | tee run_score.log
