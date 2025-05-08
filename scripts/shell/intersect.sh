#!/usr/bin/env bash

export bedtools='/data/jake/bedtools.static.binary'
gigglefile=$1
out=$2
echo "Intersecting $gigglefile with gene_file.txt.latest to $out"
"$bedtools" intersect -a /data/jake/genefusion/data/gene_file.txt.latest \
    -b <(cut -f 5- "$gigglefile") -wb -wa > "$out"


