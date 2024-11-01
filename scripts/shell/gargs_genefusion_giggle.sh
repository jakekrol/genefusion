#!/usr/bin/env bash
# parallel stix queries for pcawg prostate gene fusion data

procs=${1:-70}

#genefusion-genefusion_giggle /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt ERG 21 neg 39751949 40033704 /data/jake/genefusion/data/2024_10_31-fusions

gargs --log gargs_genefusion_giggle.log -p $procs -o "genefusion-genefusion_giggle /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt {3} {0} {4} {1} {2} /data/jake/genefusion/data/2024_10_31-fusions" < /data/jake/genefusion/data/gene_file.txt