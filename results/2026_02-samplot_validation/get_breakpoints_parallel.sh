#!/usr/bin/env bash

outdir='stix_results_summary_rnainfo_coords_giggleinfo_breakpoints'
mkdir -p "${outdir}"
tmpdir='/data/jake/tmp'
export TMPDIR=${tmpdir}
t1=$(mktemp -p "${tmpdir}" breakpoints_parallel1_XXXXXX)
trap 'rm -f "$t1"' EXIT
t2=$(mktemp -p "${tmpdir}" breakpoints_parallel2_XXXXXX)
trap 'rm -f "$t2"' EXIT
ls stix_results_summary_rnainfo_coords_giggleinfo \
    | sed 's|^|stix_results_summary_rnainfo_coords_giggleinfo/|' > $t1
ls stix_results_summary_rnainfo_coords_giggleinfo \
    | sed "s|^|${outdir}/|" > $t2

log="get_breakpoints_parallel.log"
paste $t1 $t2 | gargs -p 30 --log=$log \
    "./get_breakpoints.py -i {0} -o {1} -t ${tmpdir}"
