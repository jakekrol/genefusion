#!/usr/bin/env bash

# parallel intersect of fusion files (giggle right intervals) with gene file
# input.txt
# 0: fusionfile
# 1: outfile

export bedtools='/data/jake/bedtools.static.binary'

gargs --log gargs.log -p 60 -o "$bedtools intersect -a /data/jake/genefusion/data/gene_file.txt -b {0} > {1}" < input2.txt
