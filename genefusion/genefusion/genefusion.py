import os,sys
import time
from multiprocessing import Pool
import re
import pandas as pd
import numpy as np

def rm_non_std_chrm(df):
    allow = [str(x) for x in list(range(1,23))] + ['X', 'Y']
    mask = df.iloc[:,4].apply(lambda x: x in allow)
    return df[mask]

def rmneg1(df):
    # rm '-1' values for excord right hand hits (chr, start, end = -1); not sure what these mean
    # '-1' raises error in bedtools intersect
    filtered_df = df[(df.iloc[:, 4] != -1) & (df.iloc[:, 5] != -1) & (df.iloc[:, 6] != -1)]
    return filtered_df

def cln_giggle(df):
    df = rm_non_std_chrm(df)
    df = rmneg1(df)
    return df

def giggle_sharded(dir_shard, index, genefile, gene, chrm, strand, left, right, outdir, chdir=True):
    # giggle_sharded("/data/jake/genefusion/data/prostate/shards", "index", "/data/jake/genefusion/data/gene_file.txt", "ERG", 21, "neg", 39751949, 40033704, "/data/jake/genefusion/data/2024_10_31-fusions")
    # expect index name to be consistent across all shards
    shards = [os.path.join(dir_shard, x) for x in os.listdir(dir_shard)]
    for shard in shards:
        if chdir:
            os.chdir(shard)
        outfile = f"{os.path.basename(shard)}.{chrm}.{strand}.{gene}.{left}.{right}.giggle"
        outintersect = f"{os.path.basename(shard)}.{chrm}.{strand}.{gene}.{left}.{right}.intersect"
        outfile = os.path.join(outdir, outfile)
        outintersect = os.path.join(outdir, outintersect)
        cmd = f"/data/jake/genefusion/scripts/shell/genefusion_giggle.sh -i {index} -f {genefile} -g {gene} -c {chrm} -s {strand} -l {left} -r {right} -o {outfile} -b {outintersect}"
        os.system(cmd)

def stix(index, db, chrm, left_start, left_end, right_start, right_end, outfile,type,chdir=None):
    if chdir:
        os.chdir(chdir)
    cmd = f"stix -i {index} -d {db} -t {type} -l {chrm}:{left_start}-{left_end} -r {chrm}:{right_start}-{right_end} -s 500 > {outfile}"
    print(cmd)
    os.system(cmd)
def stix_sharded(dir_shard, chrm, left_start, left_end, right_start, right_end,outdir,outfile,type,processes,index_name='index',db_name='stix.ped.db'):
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
         data.append((index_name,db_name,chrm,left_start,left_end,right_start,right_end,shard_out,type,shard))
    print('data','\n',data[0])
    print('performing stix queries over shards at',dir_shard, 'with',processes,'processes')
    with Pool(processes=processes) as pool:
        pool.starmap(stix, data)
    print('end sharded stix')

def genefile2queries(chr_gene_file, max_dist = 5 * (10 ** 6)):
    base = os.path.basename(chr_gene_file)
    match = re.match(r"chr(\d+)\.(.+)", base)
    if match:
        chrm = match.group(1)
        strand = match.group(2)
    else:
        raise ValueError(f"file {base} does not match pattern chr.#.strand, | valid example -> chr21.neg")
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
            if gene_j_start - gene_i_stop > max_dist:
                continue
            else:
                queries['gene_i'].append(gene_i)
                queries['gene_i_start'].append(gene_i_start)
                queries['gene_i_stop'].append(gene_i_stop)
                queries['gene_j'].append(gene_j)
                queries['gene_j_start'].append(gene_j_start)
                queries['gene_j_stop'].append(gene_j_stop)

    # it's assumed these are fixed by preprocessing the input file
    queries['chr'] = [chrm] * len(queries['gene_i'])
    # strand not used for query, but use to organize output
    queries['strand'] = [strand] * len(queries['gene_i'])
    return pd.DataFrame(queries)
