#!/usr/bin/env bash

CPUS=4
mapfile -t fastqs < <(ls output/*fastq.gz)
x=$(mktemp fastq_list.XXXXXX)
trap "rm -f $x" EXIT
for f in "${fastqs[@]}"; do
	echo $f >> $x
done
seqkit stats --skip-err --seq-type rna -j $CPUS -T --infile-list $x > fastq_stats.tsv
