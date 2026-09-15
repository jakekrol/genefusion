#!/usr/bin/env python3
import numpy as np
import pandas as pd
from polymerization.score import coverage_normalize_df_evidence, coverage_normalize_tumor_reads, coverage_normalize_normal_reads, normalize_samples
import sqlite3
import time
import yaml

DATABASE="fusion.db"
TABLE='fusion_joined'
COLUMN_MAP='score_column_map.yaml'
with open(COLUMN_MAP, 'r') as f:
    column_map = yaml.safe_load(f)

### benchmark ###

# ~50k rows per second
# estimated time: ~55.5 hours total
def chunk2df(database: str, table: str, chunk_size: int = int(10**5), max_iter: int = 10):
    conn = sqlite3.connect(database)
    n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    times = []
    t_0 = time.time()
    for i,df in enumerate(pd.read_sql_query(f"SELECT * FROM {table}", conn, chunksize=chunk_size)):
        if i >= max_iter:
            break
        df_norm_cov = coverage_normalize_df_evidence(df, column_map)
        # df_norm_burden
        df_score.to_csv("chunk2df_{i}.tsv",sep="\t",index=False)
        t_elapsed = time.time() - t_0
        print(f"# time for chunk {i+1}: {t_elapsed}")
        times.append(t_elapsed)
        t_0 = time.time()
    conn.close()
    return np.mean(times), chunk_size, n


# ~125,000 rows per second
# total rows ~ 10^9
# estimated time: 22.2 hours
def chunk2py(database: str, table: str, chunk_size: int = int(10**5), max_iter: int = 10):
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
    avg_t_per_chunk, chunk_size, n  = chunk2df(DATABASE,TABLE)
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
