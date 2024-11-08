#!/usr/bin/env bash
# intersect giggle hits with gene file
set -u

bedtools='/data/jake/bedtools.static.binary'
count=false

date
while getopts "a:b:o:c" opt; do
  case $opt in
    a)
      genefile=$OPTARG
      ;;
    b)
      gigglefile=$OPTARG
      ;;
    o)
      outfile=$OPTARG
      ;;
    c)
      count=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

echo "$@"
echo "$outfile"  # This should now correctly display the outfile

if [ "$count" = "true" ]; then
    "$bedtools" intersect -a "$genefile" -b <(cut -f 5-7 "$gigglefile") | cut -f4 | sort | uniq -c | sort -n > "$outfile"
else
    "$bedtools" intersect -a "$genefile" -b <(cut -f 5-7 "$gigglefile") > "$outfile"
fi
date
