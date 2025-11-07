#!/usr/bin/env bash
set -euo pipefail

# default parameters
onekg=/data/jake/genefusion/results/2025_09-onekg_giggle2fusion/g2f/pop_normal_fusions_pe_and_sample.tsv
tumor=pop_tumor_fusions_pe_and_sample.tsv
normal=pop_normal_fusions_pe_and_sample.tsv
tumor_and_normal=pop_tumor_and_normal_fusions_pe_and_sample.tsv
all=pop_all_fusions_pe_and_sample.tsv
all_dedup=pop_all_fusions_pe_and_sample_no_dups.tsv
all_dedup_fill=pop_all_fusions_pe_and_sample_no_dups_fill.tsv

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input_dir="$2"
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
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# join tumor normal
echo "Joining tumor and normal"
join_ddb.py \
    --type outer \
    -l $tumor \
    -r $normal \
    --keys left,right \
    --output $tumor_and_normal

# join onekg
echo "Joining 1000 Genomes"
join_ddb.py -l $tumor_and_normal \
    -r  $onekg \
    --type left \
    -k left,right \
    -o $all
#de-duplicate
echo "De-duplicating"
dedup.sh \
    -i $all \
    -o $all_dedup \
    -a 1 \
    -b 2
# fillna
echo "Filling NAs with 0"
fillna.sh \
    -i $all_dedup \
    -o $all_dedup_fill \
    -f 0
