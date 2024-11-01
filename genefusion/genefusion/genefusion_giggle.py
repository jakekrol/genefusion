#!/usr/bin/env python
from genefusion.genefusion import giggle_sharded
import os,sys
# genefusion-genefusion_giggle /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt  ERG 21 neg 39751949 40033704 /data/jake/genefusion/data/2024_10_31-fusions

def main():
#    args = sys.argv[1:]
#    dir_shard = args[0]
#    index = args[1]
#    gene_file = args[2]
#    gene = args[3]
#    chrm = args[4]
#   strand = args[5]
#    left = args[5]
#    right = args[6]
#    outdir = args[7]
    file, *args = sys.argv
    print(sys.argv)
    giggle_sharded(*args)
