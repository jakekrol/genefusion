#!/usr/bin/env bash
set -euo pipefail

# default parameters
onekg=/data/jake/genefusion/results/2025_09-onekg_giggle2fusion/g2f/pop_normal_fusions_pe_and_sample.tsv
gene_pair_sorted_key=/data/jake/genefusion/results/2025_11-fusion_lexicographic_key/gene_pairs_sorted_by_genomic_position.tsv
# tumor=pop_tumor_fusions_pe_and_sample.tsv
# normal=pop_normal_fusions_pe_and_sample.tsv
# tumor_and_normal=pop_tumor_and_normal_fusions_pe_and_sample.tsv
modal=modal_combined_read_sample.tsv
all_modal=modal_onekg_combined_read_sample.tsv
all_modal_dedup=modal_onekg_combined_read_sample_no_dups.tsv
all_modal_dedup_fill=modal_onekg_combined_read_sample_no_dups_fillna.tsv

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            # a file with filenames as lines to be combined with outer join
            flist_outer="$2"
            shift 2
            ;;
        -o|--output)
            output_dir="$2"
            shift 2
            ;;
        --onekg)
            onekg="$2"
            shift 2
            ;;
        -t|--temp)
            temp_dir="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# read filenames from flist_outer
echo "Reading file list from $flist_outer"
mapfile -t files_outer < "$flist_outer"
# check each file exists
for f in "${files_outer[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "File not found: $f"
        exit 1
    fi
done

# iteratively join files in files_outer
joined_file="${temp_dir}/joined_outer.tsv"
cp "${files_outer[0]}" "$joined_file"
for ((i=1; i<${#files_outer[@]}; i++)); do
    echo "Joining file $((i+1)) of ${#files_outer[@]}"
    temp_file="${temp_dir}/temp_joined_${i}.tsv"
    join_ddb.py \
        -l "$joined_file" \
        -r "${files_outer[i]}" \
        --type outer \
        -k left,right \
        --output "$temp_file"
    mv "$temp_file" "$joined_file"
done
mv "$joined_file" "$output_dir/${modal}"

# # join tumor normal
# echo "Joining tumor and normal"
# join_ddb.py \
#     --type outer \
#     -l $tumor \
#     -r $normal \
#     --keys left,right \
#     --output $tumor_and_normal

# join onekg
echo "Joining 1000 Genomes"
join_ddb.py -l "$output_dir/$modal" \
    -r "$onekg" \
    --type left \
    -k left,right \
    -o "$output_dir/$all_modal"
#de-duplicate
echo "De-duplicating"
join_ddb.py -l $gene_pair_sorted_key \
    -r "$output_dir/$all_modal" \
    --type left \
    -k left,right \
    --output "$output_dir/$all_modal_dedup"
# fillna
echo "Filling NAs with 0"
fillna.sh \
    -i "$output_dir/$all_modal_dedup" \
    -o "$output_dir/$all_modal_dedup_fill" \
    -f 0
