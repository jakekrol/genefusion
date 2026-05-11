#!/usr/bin/env python3

import polars as pl
import time

t_0 = time.time()
df = pl.read_csv(
    'gene_pairs_sorted_by_genomic_position.tsv',
    separator='\t',
    has_header=True,
    dtypes=[pl.Utf8, pl.Utf8]
)
print(f"# tsv read time: {time.time() - t_0:.2f} seconds")
t_0=time.time()
df.write_parquet('gene_pairs_sorted_by_genomic_position.parquet')
print(f"# parquet write time: {time.time() - t_0:.2f} seconds")
