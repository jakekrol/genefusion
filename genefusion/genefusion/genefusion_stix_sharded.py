#!/usr/bin/env python
from genefusion.genefusion import stix, stix_sharded, genefile2queries
import os,sys
import multiprocessing as mp
import pandas as pd
from multiprocessing import Pool
import re
# chromosome-strand-wise genefusion stix queries

# example:  genefusion-genefusion_stix_sharded /data/jake/genefusion/data/chr_gene_lists a b c d e

def main():

    args = sys.argv[1:]
    # chromosome-strand-wise gene files
    chrgenedir = args[0]
    sharddir = args[1]
    outdir = args[2]
    outfile = args[3]
    type = args[4]
    processes = args[5]

    files = os.listdir(chrgenedir)
    files = [os.path.join(chrgenedir, f) for f in files]
    file = files[0]
    print(file)
    df = genefile2queries(file)
    print('shape', df.shape)
    print(df.head())
    print(df.tail())

        
    #def stix_sharded_genefusion(chr_gene_file,dir_shard,outdir,outfile,type,processes):
            
        
        #stix_sharded

        
        

        



