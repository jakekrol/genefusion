#!/usr/bin/env bash

# input: swapped giggle BED file
# output: intersected BED file with gene_file.txt.latest

export bedtools='/data/jake/bedtools.static.binary'
bgzipped=false

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      giggleswap="$2"
      shift 2
      ;;
    -g|--genebed)
      genebed="$2"
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

echo "Intersecting $giggleswap with $genebed to $out (bgzipped=$bgzipped)"

if [[ "$bgzipped" == true ]]; then
  # Compressed input/output: preserve comment lines and first non-comment header line
  {
    # 1) Comments + first non-comment line (column header) from swapped file
    zcat -- "$giggleswap" | awk '
      /^#/ { print; next }
      { print "gene_right_chrm\tgene_right_start\tgene_right_end\tgene_right_name\tgene_right_strand\t" $0; exit }
    '

    # 2) bedtools on data lines only (no comments, no header row in swapped file)
    "$bedtools" intersect \
      -a <(grep -v '^#' "$genebed") \
      -b <(zcat -- "$giggleswap" | awk '
          /^#/ { next }
          !header_seen { header_seen=1; next }
          { print }
        ') \
      -wb -wa
  } | bgzip -c > "$out"
else
  # Plain-text input/output: same logic without zcat/bgzip
  {
    # 1) Comments + first non-comment line (column header) from swapped file
    awk '
      /^#/ { print; next }
      { print "gene_right_chrm\tgene_right_start\tgene_right_end\tgene_right_name\tgene_right_strand\t" $0; exit }
    ' "$giggleswap"

    # 2) bedtools on data lines only
    "$bedtools" intersect \
      -a <(grep -v '^#' "$genebed") \
      -b <(awk '
          /^#/ { next }
          !header_seen { header_seen=1; next }
          { print }
        ' "$giggleswap") \
      -wb -wa
  } > "$out"
fi
