#!/usr/bin/env bash

# total
f=../../data/2025_08-haas_review_fusions/2019-haas-fusions-table-s4.txt
tail -n +2 $f | cut -f 3 | sort | uniq | wc -l > 2019-haas-fusions-table-s4-unique-gene-count.txt

# no selfies
tail -n +2 $f | cut -f 3 | sort | uniq | sed 's|--|\t|g' | awk '$1 != $2' | wc -l > 2019-haas-fusions-table-s4-unique-gene-count-no-selfies.txt