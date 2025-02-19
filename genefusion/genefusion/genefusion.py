import os, sys
import time
import multiprocessing as mp
from multiprocessing import Pool
import re
import pandas as pd
import numpy as np
from filelock import FileLock
import networkx as nx
from collections import Counter

class TreeNode:
    def __init__(self, value, index):
        self.value = value
        self.index = index  # Store index of the string in the sorted list
        self.left = None
        self.right = None

def build_balanced_tree(strings, start=0, end=None):
    if end is None:
        end = len(strings)
    if start >= end:
        return None

    mid = (start + end) // 2
    root = TreeNode(strings[mid], mid)
    root.left = build_balanced_tree(strings, start, mid)  # Left subtree
    root.right = build_balanced_tree(strings, mid + 1, end)  # Right subtree
    return root

def search_tree(root, query):
    """Search for the query string in the binary tree and return its index if found."""
    if root is None:
        return -1  # Not found
    if root.value == query:
        return root.index
    elif query < root.value:
        return search_tree(root.left, query)
    else:
        return search_tree(root.right, query)

def el2wadj_list(el, out,n = 25374):
    # n:= max gene index in /data/jake/genefusion/data/genes.index
    with open(el) as f:
        edges = [tuple(map(int, line.strip().split("\t"))) for line in f]

    # Count occurrences of each edge
    edge_counts = Counter(edges)

    # Create a weighted graph
    G = nx.Graph()
    for (u, v), weight in edge_counts.items():
        G.add_edge(u, v, weight=weight)
    # add singletons
    n = 25374
    for i in range(0, n+1):
        if i not in G.nodes:
            G.add_node(i)

    with open(out, "w") as f:
        for node, neighbors in sorted(G.adjacency(), key=lambda x: int(x[0])):
            if neighbors:
                neighbor_str = ",".join(f"({nbr}:{data['weight']})" for nbr, data in neighbors.items())
                line = f"{node},{neighbor_str}" if neighbor_str else f"{node}"
                f.write(line + "\n")
            else:
                f.write(f"{node}\n")

def index_els(el, out, index="/data/jake/genefusion/data/genes.index"):
    # input edge list of gene names
    # output edge list of gene indices
    string_list = []
    with open(index, 'r') as f:
        for line in f:
            string_list.append(line.split('\t')[1].strip('\n'))
    root = build_balanced_tree(string_list)

    with open(out, "w") as f:
        with open(el) as g:
            for line in g.readlines():
                i,j = line.strip().split("\t") # remove \n and split on tab
                g_i = i.split(".")[0] # remove version number
                g_j = j.split(".")[0] # remove version number
                # binary search to find gene index
                idx_i = search_tree(root, g_i)
                idx_j = search_tree(root, g_j)
                f.write(f"{idx_i}\t{idx_j}\n")

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
    shard_pattern="shard",
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
    shards = [os.path.join(dir_shard, x) for x in shards if re.match(shard_pattern, x)]
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


def get_sample_wise_fusions(infile, gene, outdir, append=False):
    # example
    # input:
    # infile='/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/fusions_25_01_05/10.neg.A1CF.52559169.52645435.fusion'
    # gene = 'A1CF'
    # outdir='/data/jake/genefusion/scratch/2025-01-17-sample_wise_fusions'
    # group by sample and write to file
    df = pd.read_csv(infile, sep="\t", header=None)

    # gene = os.path.basename(infile).split('.')[2:-3][0]
    print("Gene: ", gene)

    if append:
        mode = "a"
    else:
        mode = "w"

    for sample, group in df.groupby(10):  # column 11 is sample
        # strip the folder index prefix
        # and strip extension suffices
        sample = os.path.basename(sample)
        sample = sample.split(".")[0]

        group["source"] = gene
        group.rename(columns={3: "target"}, inplace=True)
        # only bother keeping source and target gene column
        group = group[["source", "target"]]
        outfile = f"{outdir}/{sample}.fusion"
        with FileLock(outfile + ".lock"):
            group.to_csv(outfile, sep="\t", header=None, index=False, mode=mode)
