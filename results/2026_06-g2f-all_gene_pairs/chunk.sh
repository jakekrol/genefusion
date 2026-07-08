#!/usr/bin/env bash

outdir="chunks"
input='all-evidence-fill.tsv'
lines_per_chunk=100000

mkdir -p "$outdir"

tail -n +2 $input | \
	split -d --lines "$lines_per_chunk" --suffix-length=4 - ${outdir}/chunk_

