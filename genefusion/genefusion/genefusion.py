import os,sys
import time
from multiprocessing import Pool
import re
import pandas as pd
import numpy as np
def stix(index, db, chr, left_start, left_end, right_start, right_end, outfile,type,chdir=None):
    if chdir:
        os.chdir(chdir)
    cmd = f"stix -i {index} -d {db} -t {type} -l {chr}:{left_start}-{left_end} -r {chr}:{right_start}-{right_end} -s 500 > {outfile}"
    print(cmd)
    os.system(cmd)
def stix_sharded(dir_shard, chr, left_start, left_end, right_start, right_end,outdir,outfile,type,processes,index_name='index',db_name='stix.ped.db'):
    print('begin sharded stix')
    processes = int(processes)
    args = locals()
    for key, value in args.items():
        print(f"{key}: {value}")
    #time.sleep(3)
    shards = os.listdir(dir_shard)
    if not shards:
        print('no shards found at',dir_shard)
        return
    shards = [os.path.join(dir_shard,shard) for shard in shards]
    data = []
    for shard in shards:
         shard_out = os.path.join(outdir,f'{os.path.basename(shard)}.{outfile}')
         data.append((index_name,db_name,chr,left_start,left_end,right_start,right_end,shard_out,type,shard))
    print('data','\n',data[0])
    print('performing stix queries over shards at',dir_shard, 'with',processes,'processes')
    with Pool(processes=processes) as pool:
        pool.starmap(stix, data)
    print('end sharded stix')

def genefile2queries(chr_gene_file):
    base = os.path.basename(chr_gene_file)
    match = re.match(r"chr(\d+)\.(.+)", base)
    if match:
        chr = match.group(1)
        strand = match.group(2)
    else:
        raise ValueError(f"file {base} does not match pattern chr(\d+)\.(.+), | valid example -> chr21.neg")
    df = pd.read_csv(chr_gene_file, sep='\t',header=None, usecols=[0,1,2,3])
    queries = {
        'gene_i': [], 'gene_i_start':[], 'gene_i_stop': [], 'gene_j': [],
        'gene_j_start': [], 'gene_j_stop': []
    }
    m = df.shape[0]
    # genefusions queries satisfy gene_i_start < gene_j_start
    for i in range(m-1):
        gene_i = df.iloc[i,3]
        gene_i_start = df.iloc[i,1]
        gene_i_stop = df.iloc[i,2]
        for k in range(i+1,m):
            gene_j = df.iloc[k,3]
            gene_j_start = df.iloc[k,1]
            gene_j_stop = df.iloc[k,2]

            queries['gene_i'].append(gene_i)
            queries['gene_i_start'].append(gene_i_start)
            queries['gene_i_stop'].append(gene_i_stop)
            queries['gene_j'].append(gene_j)
            queries['gene_j_start'].append(gene_j_start)
            queries['gene_j_stop'].append(gene_j_stop)

    queries['chr'] = [chr] * len(queries['gene_i'])
    # strand not used for query, but use to organize output
    queries['strand'] = [strand] * len(queries['gene_i'])
    return pd.DataFrame(queries)
