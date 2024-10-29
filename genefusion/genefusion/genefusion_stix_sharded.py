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
# genefusion-genefusion_stix_sharded /data/jake/genefusion/data/prostate/chr_gene_lists_split/21 /data/jake/genefusion/data/prostate/shards /data/jake/genefusion/data/prostate/fusions/21 DEL 12

def shard_match(s):
    if re.match(r"^shard", s):
        return True
    else:
        return False

def main():

    # fixed params
    args = sys.argv[1:]
    chrgenedir = args[0]
    sharddir = args[1]
    outdir = args[2]
    type = args[3]
    processes = args[4]

    if not all(list(map(os.path.exists, [chrgenedir, sharddir, outdir]))):
        print('chrgenedir', chrgenedir,'or sharddir', sharddir, 'or outdir', outdir, 'does not exist, exiting.')
        return

    files = os.listdir(chrgenedir)
    files = [os.path.join(chrgenedir, f) for f in files]
    # serially process gene files,
    # but the sharded stix queries are parallelized
    for file in files:
        df = genefile2queries(file)
        print('#query data')
        print('#shape', df.shape)
        print(df.head())
        print(df.tail())
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
            shards = [i for i in os.listdir(sharddir)]
            # pattern match filter
            shards = list(filter(shard_match, shards))
            # construct path for each shard outfile
            files_to_chk = [os.path.join(outdir, f'{i}.{outfile}') for i in shards]
            print('files_to_chk',files_to_chk)
            # check if outfile exists for all shards
            completed = list(map(os.path.exists, files_to_chk))
            if all(completed):
                print('outfile', outfile, 'already found, skipping.')
            else:
                stix_sharded(sharddir, chr, left_start, left_stop, right_start, right_stop, outdir, outfile, type, processes)      
       