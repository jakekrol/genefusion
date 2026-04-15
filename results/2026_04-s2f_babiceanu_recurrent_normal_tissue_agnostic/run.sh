#!/usr/bin/env bash

set -euo pipefail

JOIN_SCRIPT=/data/jake/rl-tools/wrangle/join_ddb.py
FINAL_TABLE=all_categories-fusion_table.tsv

# run stix2fusion pipeline
./run_s2f.py 2>&1 | tee run_s2f.log
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
