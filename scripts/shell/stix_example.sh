#!/usr/bin/env bash
# for query gene, search for gene fusions
# on same chromosome and strand
set -u
gene=${1:-TMPRSS2}
chr=${2:-21}
l_s=${3:-39751949}
l_e=${4:-40033704}
#tmp
r_s=${5:-42836478}
r_e=${6:-42903043}
echo -e "#Query\tgene: $gene\tchr: $chr" 
# need to verify fragment size for 1kg low coverage
stix -i alt_sort_b -d 1kg.ped.db -s 500 -t DEL -l $chr:$l_s-$l_e -r $chr:$r_s-$r_e
#stix -i alt_sort_b -d 1kg.ped.db -s 500 -t DEL -l 14:68603030-68603035 -r 14:68603738-68603743 | head
