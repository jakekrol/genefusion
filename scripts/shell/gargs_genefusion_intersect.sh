#!/usr/bin/env bash
# parallel intersect giggle hits with gene file
set -u

procs=${1:-64}
script="/data/jake/genefusion/scripts/shell/giggle2fusions.sh"
genefile="/data/jake/genefusion/data/gene_file.txt"

# 1kg
#input="/data/jake/genefusion/scripts/shell/giggle2fusions.input"
#outdir="/data/jake/genefusion/data/2024_11_01-fusions-1kg"
#gargs --log gargs_giggle2fusions_1kg.log -p $procs -o "$script -a $genefile -b {0} -o {1}" < $input

# pcawg
input="/data/jake/genefusion/scripts/shell/giggle2fusions-pcawg.input"
# rm -1 lines
gargs --log gargs_giggle2fusions_pcawg.log -p $procs -o "$script -a $genefile -b {0} -o {1}" < $input
