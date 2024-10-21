#!/usr/bin/env bash

dir_gene=/data/jake/gene-fusion/data/genes
# get c
mapfile -t chr_files < <(ls $dir_gene | grep "^chr")

# prefix
for i in ${!chr_files[@]}; do
    chr_files[$i]="${dir_gene}/${chr_files[$i]}"
done

echo "${chr_files[@]:0:5}"

# for each chromosome
# for each strand
# for each gene
# query stix for right-hand fusions

