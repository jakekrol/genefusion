#!/usr/bin/env bash
set -euo pipefail
t_0=$(date +%s)

script=/data/jake/genefusion/scripts/python/fids2specimen_type_counts.py
input_file=""
output_file=""

while [[ $# -gt 0 ]]; do
    case $1 in
        # file of strings that expand to each tissue's dir of beds
        # e.g., $GENEFUSION/data/dna_pad/prostate/shards/shard_*/beds
        -i|--input)
            input_file="$2"
            shift 2
            ;;
        -o|--output)
            output_file="$2"
            shift 2
            ;;
        -s|--script)
            script="$2"
            shift 2
            ;;
        -t|--tmpdir)
            export TMPDIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1
            ;;
    esac
done

# Check required arguments
if [[ -z "$input_file" ]]; then
    echo "Error: --input is required"
    exit 1
fi

if [[ -z "$output_file" ]]; then
    echo "Error: --output is required"
    exit 1
fi

# Global temp file cleanup
temp_files=()
cleanup() {
    rm -f "${temp_files[@]}"
}
trap cleanup EXIT

while IFS= read -r dir; do
    echo "# processing directory: ${dir}"
    echo "# ${dir}" >> ${output_file}
    tmp=$(mktemp specimen_counter_XXXXXX.txt)
    tmp2=$(mktemp specimen_counter2_XXXXXX.txt)
    # Add temp files to global cleanup array
    temp_files+=("$tmp" "$tmp2")
    ls ${dir} | grep -v "^$" | grep -v "^/" \
        > ${tmp}
    cut -f 1 -d '.' ${tmp} > ${tmp2} && mv ${tmp2} ${tmp}
    ${script} --fileids ${tmp} >> ${output_file}
done < "${input_file}"

t_elapsed=$(( $(date +%s) - t_0 ))
echo "# done. Elapsed time: ${t_elapsed} seconds."
echo "# output written to: ${output_file}"