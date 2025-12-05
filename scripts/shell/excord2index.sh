#!/usr/bin/env bash

set -euo pipefail

t_0=$(date +%s)

SCRIPT_SHARD="/data/jake/genefusion/scripts/python/giggle_shard.py"
SCRIPT_CLN="/data/jake/genefusion/scripts/shell/cln_excord.sh"
SCRIPT_SORT="/data/jake/genefusion/scripts/shell/sort_bed"
SCRIPT_SINGULARITY_ENTRY="$HOME/sv/run.sh"

# args
echo "# parsing arguments"
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input_beds="$2"
            shift 2
            ;;
        -o|--output)
            output_dir="$2"
            shift 2
            ;;
        -n|--name)
            name="$2"
            shift 2
            ;;
        -k|--num_shards)
            num_shards="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

echo "# verifying required arguments"
for arg in input_beds output_dir name num_shards; do
    if [[ -z "${!arg:-}" ]]; then
        echo "Error: Missing required argument: $arg"
        exit 1
    fi
done

# make output dirs
echo "# making output directories"
mkdir -p "$output_dir"
mkdir -p "$output_dir/beds"
mkdir -p "$output_dir/beds_cln"
mkdir -p "$output_dir/beds_sort_cln"
mkdir -p "$output_dir/shards"

# copy input beds to output_dir/beds
echo "# copying input beds to output directory"
cp "$input_beds"/* "$output_dir/beds"

# clean beds
echo "# cleaning beds"
cd "$output_dir/beds" || { echo "Failed to cd to $output_dir/beds"; exit 1; }
ls | gargs --log "$output_dir/clean_excord.log" \
    -p 60 \
    -o "$SCRIPT_CLN -i {0} -o $output_dir/beds_cln/{0} -z"

# sort beds
echo "# sorting beds"
cd "$output_dir" || { echo "Failed to cd to $output_dir"; exit 1; }
$SCRIPT_SORT "$output_dir/beds_cln" "$output_dir/beds_sort_cln" 50 2>&1 | tee "$output_dir/sort.log"

# shard beds
echo "# sharding beds"
$SCRIPT_SHARD \
        -i "$output_dir/beds_sort_cln" \
        -o "$output_dir/shards" \
        -n "$num_shards" -s "$output_dir/shards.yaml" 2>&1 | tee "$output_dir/shard.log"

# enter singularity
echo "# entering singularity"
$SCRIPT_SINGULARITY_ENTRY

echo "# creating sharded indices"
# go back to output dir
cd "$output_dir" || { echo "Failed to cd to $output_dir"; exit 1; }
# index shards
cd shards || { echo "Failed to cd to $output_dir/shards"; exit 1; }
# beds is a dir within each shard
input='"beds/*.gz"'
# make indices in parallel
seq 0 $((num_shards - 1)) | sed "s|^|shard_|" | \
    gargs -p "$num_shards" --log=index.log -o \
    "cd {0} && giggle index -s -i $input -o ${name}_sort_b"

t_end=$(date +%s)
duration=$((t_end - t_0))
echo "# completed in $duration seconds"