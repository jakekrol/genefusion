#!/usr/bin/env bash

# input: intersect BED file with 
# - c1-5 is from genefile
#         - c6-9 is rh
#         -  c10-13 is lh
# - c14 is ?
# - c15 is sample
# output: intersect BED file with lh at c6-9 and rh at c10-13

bgzipped=false

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

echo "Unswapping intervals from $in to $out (bgzipped=$bgzipped)"

if [[ "$bgzipped" == true ]]; then
  zcat -- "$in" \
  | awk 'BEGIN {OFS="\t"}
         /^#/ {print; next}
         {print $1, $2, $3, $4, $5, $10, $11, $12, $13, $6, $7, $8, $9, $14, $15}' \
  | bgzip -c > "$out"
else
  awk 'BEGIN {OFS="\t"}
       /^#/ {print; next}
       {print $1, $2, $3, $4, $5, $10, $11, $12, $13, $6, $7, $8, $9, $14, $15}' \
      "$in" > "$out"
fi
