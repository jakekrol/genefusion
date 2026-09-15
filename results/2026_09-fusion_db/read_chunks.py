#!/usr/bin/env python3
import numpy as np
import pandas as pd
from polymerization.score import coverage_normalize_df_evidence, coverage_normalize_tumor_reads, coverage_normalize_normal_reads, normalize_samples, burden_normalize_df_evidence, burden_normalize_reads
import sqlite3
import time
import yaml

DATABASE="fusion.db"
TABLE='fusion_joined'
TABLE_BURDEN='burden'
COLUMN_MAP='score_column_map.yaml'
with open(COLUMN_MAP, 'r') as f:
    column_map = yaml.safe_load(f)

### benchmark ###

def chunk2df(database: str, table: str, table_burden: str, chunk_size: int = int(10**5), max_iter: int = 10):
    # burden
    conn = sqlite3.connect(database)
    df_burden = pd.read_sql_query(f"SELECT * FROM {table_burden}", conn)
    df_burden.set_index("gene", inplace=True)
    # evidence
    n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    times = []
    t_0 = time.time()
    for i,df in enumerate(pd.read_sql_query(f"SELECT * FROM {table}", conn, chunksize=chunk_size)):
        if i >= max_iter:
            break
        df_norm_coverage = coverage_normalize_df_evidence(df.copy(), column_map)
        df_norm_burden = burden_normalize_df_evidence(df, df_burden, column_map)
        df_norm_coverage.to_csv(f"chunk2df_norm_cov_{i}.tsv",sep="\t",index=False)
        df_norm_burden.to_csv(f"chunk2df_norm_burden_{i}.tsv",sep="\t",index=False)
        t_elapsed = time.time() - t_0
        print(f"# time for chunk {i+1}: {t_elapsed}")
        times.append(t_elapsed)
        t_0 = time.time()
    conn.close()
    return np.mean(times), chunk_size, n


def chunk2py(database: str, table: str, table_burden: str, chunk_size: int = int(10**5), max_iter: int = 10):
    conn = sqlite3.connect(database)
    n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    times = []
    cursor = conn.execute(f"SELECT * FROM {table}")
    i=1
    t_0 = time.time()
    while (rows := cursor.fetchmany(chunk_size)) and (i < 11):
        row_count = i * chunk_size
        for row in rows:
            len(row)
        print(f"row {row_count}/{n}")
        t_elapsed = time.time() - t_0
        print(f"# time for chunk {i+1}: {t_elapsed}")
        times.append(t_elapsed)
        t_0 = time.time()
        i+=1
    conn.close()
    return np.mean(times), chunk_size, n

def estimate_runtime(avg_t_per_chunk, chunk_size, n):
    iterations = n / chunk_size
    t_est = avg_t_per_chunk * iterations
    return t_est

def benchmark():
    # estimated time: ~55.5 hours total
    print("# benchmarking chunk2df")
    avg_t_per_chunk, chunk_size, n  = chunk2df(DATABASE,TABLE,TABLE_BURDEN)
    t_est = estimate_runtime(avg_t_per_chunk, chunk_size, n)
    print(
        "# estimated time chunk2df with rows={} and chunksize={}: {}".format(
            n, chunk_size, t_est
        )
    )
    # estimated time: 22.2 hours
    print("# benchmarking chunk2py")
    avg_t_per_chunk, chunk_size, n  = chunk2py(DATABASE,TABLE)
    t_est = estimate_runtime(avg_t_per_chunk, chunk_size,n)
    print(
        "# estimated time chunk2py with rows={} and chunksize={}: {} seconds".format(
            n, chunk_size, t_est
        )
    )


def main():
    # results: chunk2py=3356 seconds and chunk2df=9805 seconds
    benchmark()

if __name__ == "__main__":
    main()
