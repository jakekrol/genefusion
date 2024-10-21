#!/usr/bin/env python
from genefusion.genefusion import stix, stix_sharded, genefile2queries
import os,sys
import multiprocessing as mp
import pandas as pd
from multiprocessing import Pool
import re
# chromosome-strand-wise genefusion stix queries

# example:  
# genefusion-genefusion_stix_sharded /data/jake/genefusion/data/chr_gene_lists /data/jake/genefusion/data/prostate/shards /data/jake/genefusion/data/fusions_test DEL 20

def main():

    # fixed params
    args = sys.argv[1:]
    chrgenedir = args[0]
    sharddir = args[1]
    outdir = args[2]
    type = args[3]
    processes = args[4]

    files = os.listdir(chrgenedir)
    files = [os.path.join(chrgenedir, f) for f in files]
    # serially process gene files,
    # but the sharded stix queries are parallelized
    for file in files:
        df = genefile2queries(file)
        print('shape', df.shape)
        print(df.head())
        print(df.tail())
        # genefusion query params
        for idx, row in df.iterrows():
            chr = row['chr']
            strand = row['strand']
            gene_i = row['gene_i']
            left_start = row['gene_i_start']
            left_stop = row['gene_i_stop']
            right_start = row['gene_j_start']
            right_stop = row['gene_j_stop']
            gene_j = row['gene_j']
            outfile = f'{chr}.{strand}.{gene_i}.{gene_j}.stix_sharded.out'
            stix_sharded(sharddir, chr, left_start, left_stop, right_start, right_stop, outdir, outfile, type, processes)
        