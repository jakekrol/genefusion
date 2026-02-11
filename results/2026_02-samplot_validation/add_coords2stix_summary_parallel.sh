#!/usr/bin/env bash

din=stix_results_summary
dout=stix_results_summary_rnainfo_coords
logfile=add_coords2stix_summary_parallel.log
procs=30
mkdir -p $dout
t1=$(mktemp add_coords2stix_summary_parallel.XXXXXX)
trap 'rm -f "$t1"' EXIT
t2=$(mktemp add_coords2stix_summary_parallel2.XXXXXX)
trap 'rm -f "$t2"' EXIT
ls $din/*.out > $t1
ls $din > $t2
sed -i "s|^|$dout/|" $t2
paste $t1 $t2 | \
    gargs --log=$logfile -p $procs \
        "LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH ./add_coords2stix_summary.py -i {0} -o {1}"

