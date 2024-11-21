#!/usr/bin/env bash

# parallel intersect of (giggle right intervals) with gene file
# input.txt
# 0: gigglefile (right intervals cut)
# 1: outfile

export bedtools='/data/jake/bedtools.static.binary'
# ls | gargs -o -p 60 "cut -f 5-7 {0} > ../giggle_cut/{0}"
# for f in $(ls giggle_cut); do base=$(echo ${f%.*}); printf "giggle_cut/$f\t fusions2/${base}.fusions\n" >> input.txt; done

gargs --log gargs.log -p 64 -o "$bedtools intersect -a /data/jake/genefusion/data/gene_file.txt -b {0} > {1}" < input.txt
