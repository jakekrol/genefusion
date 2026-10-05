#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import sqlite3
import subprocess
from multiprocessing import Pool
from polymerization.datasets import get_pcawg_data_types

parser=argparse.ArgumentParser()
parser.add_argument("--db", default="fusion.db")
parser.add_argument("--outdir", default="top_fusions")
parser.add_argument("--table", default="score")
parser.add_argument("--num_fusions", default=10000)
parser.add_argument("--cpus",default=30)
args = parser.parse_args()
args.num_fusions = int(args.num_fusions)
args.cpus = int(args.cpus)
os.makedirs(args.outdir,exist_ok=True)

def get_pcawg_tumor_tissues():
    data = get_pcawg_data_types()
    tissues = set()
    for tissue, _ in data.items():
        if 'tumor' in data[tissue]['modality']['dna']['specimen']:
            tissues.add(tissue)
        if 'tumor' in data[tissue]['modality']['rna']['specimen']:
            tissues.add(tissue)
    return tissues

def write_topk(db, table, column, outfile, num_fusions=1000):
    conn = sqlite3.connect(db)
    sql = f"SELECT gene_left, gene_right, {column} FROM {table} ORDER BY {column} DESC LIMIT {num_fusions};"
    print(f"# sql query: {sql}")
    df = pd.read_sql_query(sql, conn)
    df.to_csv(outfile,sep="\t", index=False)
    return df

    

def main():
    df_global_norm_coverage = write_topk(
        args.db,
        args.table,
        "score_global_norm_coverage",
        os.path.join(args.outdir,f"top-{args.num_fusions}-score-global_norm_coverage.tsv"),
        args.num_fusions
    )
    df_global_norm_burden = write_topk(
        args.db,
        args.table,
        "score_global_norm_burden",
        os.path.join(args.outdir,f"top-{args.num_fusions}-score-global_norm_burden.tsv"),
        args.num_fusions
    )
    queries = []
    tumor_tissues = get_pcawg_tumor_tissues()
    for tissue in tumor_tissues:
        for score_type in ["coverage", "burden"]:
            column = f"score_{tissue}_norm_{score_type}"
            outfile = f"top-{args.num_fusions}-score-{tissue}-{score_type}.tsv"
            outfile = os.path.join(
                args.outdir,
                outfile
            )
            queries.append((args.db, args.table, column, outfile, args.num_fusions))
    n = len(queries)
    print(f"# total queries: {n}")
    cpus = min(args.cpus, n)
    with Pool(cpus) as pool:
        results = pool.starmap(write_topk, queries)

if __name__ == "__main__":
    main()