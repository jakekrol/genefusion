#!/usr/bin/env bash

dirinput='./stix_results_summary'
diroutput='./stix_results_summary_rnainfo'
mkdir -p $diroutput
t1=$(mktemp "add_paired_rna_filename_parallel_XXXXXX.tsv")
trap 'rm -f $t1' EXIT
t2=$(mktemp "add_paired_rna_filename_parallel2_XXXXXX.tsv")
trap 'rm -f $t2' EXIT
ls $dirinput | sed "s|^|$dirinput/|" > $t1
ls $dirinput | sed "s|^|$diroutput/|" > $t2
paste $t1 $t2 | \
    gargs --log="add_paired_rna_filename_parallel.log" \
        -p 30 \
        "./add_paired_rna_filename.py -i {0} -o {1}"
