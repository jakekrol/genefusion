import os, sys
import time
import multiprocessing as mp
from multiprocessing import Pool
import re
import pandas as pd
import numpy as np


def rm_non_std_chrm(df):
    allow = [str(x) for x in list(range(1, 23))] + ["X", "Y"]
    mask = df.iloc[:, 4].apply(lambda x: x in allow)
    return df[mask]


def rmneg1(df):
    # rm '-1' values for excord right hand hits (chr, start, end = -1); not sure what these mean
    # '-1' raises error in bedtools intersect
    filtered_df = df[
        (df.iloc[:, 4] != -1) & (df.iloc[:, 5] != -1) & (df.iloc[:, 6] != -1)
    ]
    return filtered_df


def rm_double_zero(df):
    # rm rows where both start and end are 0
    rm_idx = df[(df.iloc[:, 5] == 0) & (df.iloc[:, 6] == 0)].index
    if len(rm_idx) > 0:
        return df.drop(rm_idx)
    else:
        return df


def cln_giggle(df):
    """
    clean giggle files by removing non-standard chromosomes, -1 values for chromosomes,
    and removes rows where start and end are both 0
    from giggle output
    """
    df = rm_non_std_chrm(df)
    df = rmneg1(df)
    df = rm_double_zero(df)
    return df


def giggle_sharded(
    dir_shard,
    index,
    genefile,
    gene,
    chrm,
    strand,
    left,
    right,
    outdir,
    chdir=True,
    # intersect=False,
    shard_pattern="shard",
    parallel=False,
    cpus=mp.cpu_count(),
):
    """
    dir_shard: str, full path to directory containing shards
    index: str, name of giggle index should be the same for each shard
    genefile: str, full path to gene file
    gene: str, gene name
    chrm: str, chromosome
    strand: str, strand
    left: int, left coordinate
    right: int, right coordinate
    outdir: str, full path to output directory
    chdir: bool, change directory to shard before running giggle
    shard_pattern: str, pattern to match shards
    parallel: bool, parallelize giggle calls
    cpus: int, number of cpus to use
    """
    cpus = int(cpus)

    # setup paths
    shards = os.listdir(dir_shard)
    shards = [os.path.join(dir_shard, x) for x in shards if re.match(shard_pattern, x)]
    outfiles = [
        f"{os.path.basename(x)}.{chrm}.{strand}.{gene}.{left}.{right}.giggle"
        for x in shards
    ]
    outfiles = [os.path.join(outdir, x) for x in outfiles]
    # setup calls to shell script
    cmds = [
        f"/data/jake/genefusion/scripts/shell/genefusion_giggle.sh -i \
            {index} -f {genefile} -g {gene} -c {chrm} -s {strand} -l {left} -r {right} -o {outfile}"
        for outfile in outfiles
    ]
    # if intersect:
    #     cmds = [
    #         f"{cmd} -b {outfile.replace('.giggle', '.intersect')}"
    #         for cmd, outfile in zip(cmds, outfiles)
    #     ]
    # optionally parallelize
    if not parallel:
        for cmd, shard in zip(cmds, shards):
            if chdir:
                os.chdir(shard)
            os.system(cmd)
    else:
        cpus = min(cpus, len(shards))
        cmds = [f"cd {shard}; {cmd}" for cmd, shard in zip(cmds, shards)]
        [print(cmd) for cmd in cmds]
        with Pool(processes=cpus) as pool:
            pool.map(os.system, cmds)

def stix(
    index,
    db,
    chrm,
    left_start,
    left_end,
    right_start,
    right_end,
    outfile,
    type,
    chdir=None,
):
    if chdir:
        os.chdir(chdir)
    cmd = f"stix -i {index} -d {db} -t {type} -l {chrm}:{left_start}-{left_end} -r {chrm}:{right_start}-{right_end} -s 500 > {outfile}"
    print(cmd)
    os.system(cmd)


def stix_sharded(
    dir_shard,
    chrm,
    left_start,
    left_end,
    right_start,
    right_end,
    outdir,
    outfile,
    type,
    processes,
    index_name="index",
    db_name="stix.ped.db",
):
    print("begin sharded stix")
    processes = int(processes)
    args = locals()
    for key, value in args.items():
        print(f"{key}: {value}")
    # time.sleep(3)
    shards = os.listdir(dir_shard)
    if not shards:
        print("no shards found at", dir_shard)
        return
    shards = [os.path.join(dir_shard, shard) for shard in shards]
    data = []
    for shard in shards:
        shard_out = os.path.join(outdir, f"{os.path.basename(shard)}.{outfile}")
        data.append(
            (
                index_name,
                db_name,
                chrm,
                left_start,
                left_end,
                right_start,
                right_end,
                shard_out,
                type,
                shard,
            )
        )
    print("data", "\n", data[0])
    print(
        "performing stix queries over shards at",
        dir_shard,
        "with",
        processes,
        "processes",
    )
    with Pool(processes=processes) as pool:
        pool.starmap(stix, data)
    print("end sharded stix")


def genefile2queries(chr_gene_file, max_dist=5 * (10**6)):
    base = os.path.basename(chr_gene_file)
    match = re.match(r"chr(\d+)\.(.+)", base)
    if match:
        chrm = match.group(1)
        strand = match.group(2)
    else:
        raise ValueError(
            f"file {base} does not match pattern chr.#.strand, | valid example -> chr21.neg"
        )
    df = pd.read_csv(chr_gene_file, sep="\t", header=None, usecols=[0, 1, 2, 3])
    queries = {
        "gene_i": [],
        "gene_i_start": [],
        "gene_i_stop": [],
        "gene_j": [],
        "gene_j_start": [],
        "gene_j_stop": [],
    }
    m = df.shape[0]
    # genefusions queries satisfy gene_i_start < gene_j_start
    for i in range(m - 1):
        gene_i = df.iloc[i, 3]
        gene_i_start = df.iloc[i, 1]
        gene_i_stop = df.iloc[i, 2]
        for k in range(i + 1, m):
            gene_j = df.iloc[k, 3]
            gene_j_start = df.iloc[k, 1]
            gene_j_stop = df.iloc[k, 2]
            if gene_j_start - gene_i_stop > max_dist:
                continue
            else:
                queries["gene_i"].append(gene_i)
                queries["gene_i_start"].append(gene_i_start)
                queries["gene_i_stop"].append(gene_i_stop)
                queries["gene_j"].append(gene_j)
                queries["gene_j_start"].append(gene_j_start)
                queries["gene_j_stop"].append(gene_j_stop)

    # it's assumed these are fixed by preprocessing the input file
    queries["chr"] = [chrm] * len(queries["gene_i"])
    # strand not used for query, but use to organize output
    queries["strand"] = [strand] * len(queries["gene_i"])
    return pd.DataFrame(queries)
