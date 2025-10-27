#!/usr/bin/env bash

awk -v OFS="\t" '{print $0, $4"."$1"."$5"."$2"."$3}' \
    ../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed \
    > giggle_search.template
