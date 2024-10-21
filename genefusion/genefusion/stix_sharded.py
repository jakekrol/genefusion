#!/usr/bin/env python
from genefusion.genefusion import stix, stix_sharded
import os,sys
import multiprocessing as mp
from multiprocessing import Pool
# example
# genefusion-stix_sharded /data/jake/genefusion/data/prostate/shards 21 39751949 40033704 42836478 42903043 /data/jake/genefusion/data/fusions erg.tmprss2.stix_sharded.out DEL 8

def main():
    args = sys.argv[1:]
    dir_shard = args[0]
    chr = args[1]
    left_start = args[2]
    left_end = args[3]
    right_start = args[4]
    right_end = args[5]
    outdir = args[6]
    outfile = args[7]
    type = args[8]
    proccesses = args[9]

    file,*stix_sharded_args = sys.argv

    print(stix_sharded_args)
    stix_sharded(*stix_sharded_args)



