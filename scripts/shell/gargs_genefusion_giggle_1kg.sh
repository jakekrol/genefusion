#!/usr/bin/env bash
# parallel stix queries for pcawg prostate gene fusion data

procs=${1:-70}
1kg="/data/jake/genefusion/data/1kg-lc"
cd $1kg || exit 1
script="/data/jake/genefusion/scripts/shell/genefusion_giggle.sh"
index='alt_sort_b'
#index="/data/jake/genefusion/data/1kg-lc/alt_sort_b"
genefile="/data/jake/genefusion/data/gene_file.txt"
outdir="/data/jake/genefusion/data/2024_11_01-fusions-1kg"
#genefusion-genefusion_giggle /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt ERG 21 neg 39751949 40033704 /data/jake/genefusion/data/2024_10_31-fusions

gargs --log gargs_genefusion_giggle_1kg.log -p $procs -o "$script -i $index -f /data/jake/genefusion/data/gene_file.txt -c {0} -s {4} -g {3} -l {1} -r {2} -o $outdir/{0}.{4}.{3}.{1}.{2}.giggle -b $outdir/{0}.{4}.{3}.{1}.{2}.intersect < /data/jake/genefusion/data/gene_file.txt