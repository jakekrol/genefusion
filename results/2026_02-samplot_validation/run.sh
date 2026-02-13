#!/usr/bin/env bash

# get top candidates from score files
./get_top_candidates.sh

# setup stix queries
for tissue in kidney liver ovary; do
    input_file="./${tissue}_top_100_candidates.tsv"
    output_file="./${tissue}_top_candidates-stix_coords.tsv"
    echo "# processing tissue: ${tissue}"
    echo "# input file: ${input_file}"
    echo "# output file: ${output_file}"
    ./fusion2stix_query.sh -i <(cut -f 1,2 "${input_file}") -o "${output_file}"
done

# run stix queries
# OUTDIR="$(pwd)/stix_results"
./run_stix.sh


# total the read evidence for each sample and get top samples for each fusion
# diroutput='./stix_results_summary'
./summarize_stix_results_parallel.sh

# add rna file id mapping here
# and filter out samples without corresponding RNA data
# diroutput='./stix_results_summary_rnainfo'
./add_paired_rna_filename_parallel.sh

# add coordinates to the top candidates to run giggle
# need to run giggle for breakpoints
# dout=stix_results_summary_rnainfo_coords
./add_coords2stix_summary_parallel.sh

# add more info needed for giggle queries
# outdir is stix_results_sum_coords_giggleinfo
mkdir -p stix_results_summary_rnainfo_coords_giggleinfo
ls stix_results_summary_rnainfo_coords \
    | sed 's|^|stix_results_summary_rnainfo_coords/|' > x
ls stix_results_summary_rnainfo_coords \
    | sed 's|^|stix_results_summary_rnainfo_coords_giggleinfo/|' > y
log="stix_summary2giggleinfo.log"
paste x y | gargs -p 30 --log=$log \
    "./stix_summary2giggleinfo.py -i {0} -o {1}"
rm x y

# run giggle queries and get breakpoints
# !note: need to have giggle in path
./get_breakpoints_parallel.sh

# get download regions for each breakpoint
./set_download_region_parallel.sh

# combine all plotting data into single table
./combine_all_plot_data.py -i stix_results_summary_rnainfo_coords_giggleinfo_breakpoints_region \
    -o final_plot_data.tsv


