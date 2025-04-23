#!/usr/bin/env bash
# clean excord files
in=$1
out=$2
echo "cleaning $in"

awk '$1 != "*"' $in | \
    awk '$2 != "0" && $3 != "0"' | \
    awk '$6 != "0" && $7 != "0"' | \
    awk '$1 !~ /^hs/ && $5 !~ /^hs/' | \
    awk '$1 !~ /^GL/ && $5 !~ /^GL/' | \
    awk '$1 !~ /^NC/ && $5 !~ /^NC/' | \
    awk '$1 !~ /^MT/ && $5 !~ /^MT/' | \
    awk '$1 !~ "-1" && $5 !~ "-1"' > "$out"
    # bgzip -c > "$out"
