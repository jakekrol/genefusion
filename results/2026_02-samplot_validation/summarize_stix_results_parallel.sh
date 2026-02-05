#!/usr/bin/env bash

k=3
dirinput='./stix_results'
diroutput='./stix_results_summary'
mkdir -p $diroutput
t1=$(mktemp "summarize_stix_results_parallel_XXXXXX.tsv")
trap 'rm -f $t1' EXIT
t2=$(mktemp "summarize_stix_results_parallel2_XXXXXX.tsv")
trap 'rm -f $t2' EXIT
ls $dirinput | sed "s|^|$dirinput/|" > $t1
ls $dirinput | sed "s|^|$diroutput/top_${k}.|" > $t2
paste $t1 $t2 | \
    gargs --log="summarize_stix_results_parallel.log" \
        -p 30 \
        "./summarize_stix_results.py -i {0} -o {1} -k ${k}"
