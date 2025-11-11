#!/usr/bin/env bash
# clean excord files
while getopts "i:o:z" opt; do
  case $opt in
    i) in="$OPTARG"
    ;;
    o) out="$OPTARG"
    ;;
    # is input bgzipped?
    z) bgzipped=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
echo "cleaning $in"

if [ -z "$in" ] || [ -z "$out" ]; then
    echo "Usage: $0 -i input_excord.tsv -o output_cleaned_excord.tsv"
    exit 1
fi

# if bgzipped is true, then use zcat and bgzip the output
if [ "$bgzipped" = true ]; then
    zcat "$in" | awk '$1 != "*"' | \
        awk '$2 != "0" && $3 != "0"' | \
        awk '$6 != "0" && $7 != "0"' | \
        awk '$1 !~ /^hs/ && $5 !~ /^hs/' | \
        awk '$1 !~ /^GL/ && $5 !~ /^GL/' | \
        awk '$1 !~ /^NC/ && $5 !~ /^NC/' | \
        awk '$1 !~ /^MT/ && $5 !~ /^MT/' | \
        awk '$1 !~ "-1" && $5 !~ "-1"' | \
        bgzip -c > "$out"
else 
    awk '$1 != "*"' "$in" | \
        awk '$2 != "0" && $3 != "0"' | \
        awk '$6 != "0" && $7 != "0"' | \
        awk '$1 !~ /^hs/ && $5 !~ /^hs/' | \
        awk '$1 !~ /^GL/ && $5 !~ /^GL/' | \
        awk '$1 !~ /^NC/ && $5 !~ /^NC/' | \
        awk '$1 !~ /^MT/ && $5 !~ /^MT/' | \
        awk '$1 !~ "-1" && $5 !~ "-1"' > "$out"
fi