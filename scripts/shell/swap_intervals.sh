#!/usr/bin/env bash

# input: giggle BED file
# output: BED file with left and right intervals swapped (c1-4 -> c5-8 & c5-8 -> c1-4)
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      in="$2"
      shift 2
      ;;
    -o|--output)
      out="$2"
      shift 2
      ;;
    # is input bgzipped?
    -z|--bgzipped)
      bgzipped=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

if [[ $bgzipped == true ]]; then
    echo "Swapping intervals in $in to $out (bgzipped)"
    zcat -- "$in" \
    | awk 'BEGIN {OFS="\t"}
           /^#/ {print; next}
           { print $5, $6, $7, $8, $1, $2, $3, $4, $9, $10 }' \
    | bgzip -c > "$out"
else
    echo "Swapping intervals in $in to $out"
    awk 'BEGIN {OFS="\t"}
         /^#/ {print; next}
         { print $5, $6, $7, $8, $1, $2, $3, $4, $9, $10 }' \
        "$in" > "$out"
fi