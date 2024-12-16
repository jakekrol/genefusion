#!/usr/bin/env python
from genefusion.genefusion import giggle_sharded
import os, sys

### new
# serial
# gf-giggle_sharded /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt  ERG 21 neg 39751949 40033704 /data/jake/genefusion/scratch/2024-12-15-giggle_sharded
# parallel
# gf-giggle_sharded /data/jake/genefusion/data/prostate/shards index /data/jake/genefusion/data/gene_file.txt  ERG 21 neg 39751949 40033704 /data/jake/genefusion/scratch/2024-12-15-giggle_sharded True shard True 4

def main():
    file, *args = sys.argv
    print(sys.argv)
    giggle_sharded(*args)
