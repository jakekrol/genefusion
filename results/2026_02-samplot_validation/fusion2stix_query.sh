#!/usr/bin/env bash

# receive as input a 2 column file of fusion gene pairs (tab delimited)
# lookup the coordinates of each gene from a gene annotation file
t_0=$(date +%s)
bedfile="../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed"
tmpdir=/data/jake/tmp
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input_file="$2"
            shift 2
            ;;
        -o|--output)
            output_file="$2"
            shift 2
            ;;
        -b|--bedfile)
            bedfile="${2:-../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed}"
            shift 2
            ;;
        -t|--tempdir)
            temp_dir="$2"
            export TMPDIR="$temp_dir"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done


if [[ -z "${input_file}" || -z "${output_file}" ]]; then
    echo "Usage: $0 -i <input_file> -o <output_file> [-f] [-b <bedfile>]"
    exit 1
fi

if [[ ! -e "${input_file}" && ! -r "${input_file}" ]]; then
    echo "Error: input file not found or not readable: ${input_file}"
    exit 1
fi
if [[ ! -f "${bedfile}" ]]; then
    echo "Error: bedfile not found: ${bedfile}"
    exit 1
fi

echo "# begin fusion2stix_query.sh for input: ${input_file}, output: ${output_file}, bedfile: ${bedfile}, temp_dir: ${temp_dir}"
tmp_out=$(mktemp "${TMPDIR:-/tmp}/fusion2stix_query_XXXXXX.tsv")
trap 'rm -f "${tmp_out}"' EXIT
while IFS=$'\t' read -r gene1 gene2; do
    coord1=$(grep -w "${gene1}" "${bedfile}" | head -n 1 | cut -f 1-3)
    coord2=$(grep -w "${gene2}" "${bedfile}" | head -n 1 | cut -f 1-3)
    if [[ -n "${coord1}" && -n "${coord2}" ]]; then
        end1=$(echo "${coord1}" | awk -v FS="\t" '{print $3}')
        start2=$(echo "${coord2}" | awk -v FS="\t" '{print $2}')
        if [[ "${end1}" =~ ^[0-9]+$ && "${start2}" =~ ^[0-9]+$ ]]; then
            # because gene left (gene1) is always smaller position we can compute distance
            dist=$((start2 - end1))
        else
            dist=""
        fi
        printf "%s\t%s\t%s\t%s\t%s\n" "${gene1}" "${coord1}" "${gene2}" "${coord2}" "${dist}"
    else
        echo "Warning: Gene not found in bedfile: ${gene1} or ${gene2}" >&2
    fi
done < "${input_file}" > "${tmp_out}"

# filter out overlapping genes by checking if distance is negative
awk '{
    if ($9 >= 0) {
        print $0
    }
}' "${tmp_out}" > "${output_file}"

t_1=$(date +%s)
echo "# completed in $((t_1 - t_0)) seconds. output written to: ${output_file}"
