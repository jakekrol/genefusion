#!/usr/bin/env bash

outdir='stix_results_summary_rnainfo_coords_giggleinfo_breakpoints_region'
mkdir -p "${outdir}"
tmpdir='/data/jake/tmp'
export TMPDIR=${tmpdir}
t1=$(mktemp -p "${tmpdir}" region_parallel1_XXXXXX)
trap 'rm -f "$t1"' EXIT
t2=$(mktemp -p "${tmpdir}" region_parallel2_XXXXXX)
trap 'rm -f "$t2"' EXIT
ls stix_results_summary_rnainfo_coords_giggleinfo_breakpoints \
    | sed 's|^|stix_results_summary_rnainfo_coords_giggleinfo_breakpoints/|' > $t1
ls stix_results_summary_rnainfo_coords_giggleinfo_breakpoints \
    | sed "s|^|${outdir}/|" > $t2

log="set_download_region_parallel.log"
paste $t1 $t2 | gargs -p 30 --log=$log \
    "./set_download_region.py -i {0} -o {1} --len_multiplier 2.0"