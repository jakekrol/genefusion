#!/usr/bin/env python
from genefusion.genefusion import stix, stix_sharded
import os,sys
import multiprocessing as mp
from multiprocessing import Pool
# example
# genefusion-stix_sharded /data/jake/genefusion/data/prostate/shards 21 39751949 40033704 42836478 42903043 /data/jake/genefusion/data/fusions erg.tmprss2.stix_sharded.out DEL 8

def main():
    file,*stix_sharded_args = sys.argv
    print(stix_sharded_args)
    stix_sharded(*stix_sharded_args)



