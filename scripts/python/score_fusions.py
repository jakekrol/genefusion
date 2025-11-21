#!/usr/bin/env python3
import duckdb
# import swifter
# import dask
import argparse
import numpy as np
from genefusion.genefusion import score_numba_vectorized
import yaml
import os,sys
import time
import numba

parser = argparse.ArgumentParser(description='Score fusions')
parser.add_argument('-i', '--input', type=str,required=True, help='Input fusion feature file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output scored fusion file')
parser.add_argument('--score_yaml',type=str,help='YAML file with population sizes and column names', required=True)
# parser.add_argument('--batch_size',type=int,default=20000000,help='Number of rows to process per batch (default: 10000000)')
parser.add_argument('--batch_size',type=int,default=80000000,help='Number of rows to process per batch')
parser.add_argument('--n_threads',type=int,default=os.cpu_count()-10,help='Number of CPUs to use (default: all available)')
args = parser.parse_args()

conn = duckdb.connect()
total_rows = conn.execute(f'SELECT COUNT(*) FROM read_csv_auto("{args.input}", delim="\t")').fetchone()[0]
print(f"Total rows: {total_rows}")


with open(args.score_yaml) as f:
    params = yaml.safe_load(f)
# print out params for debugging
print("Scoring parameters:")
for k,v in params.items():
    print(f"  {k}: {v}")

# Set numba thread count
numba.set_num_threads(args.n_threads)
print(f"Numba configured to use {args.n_threads} threads")

# warm up numba with small array to compile and test parallel paths
print("Warming up numba...")
warmup_size = 10
_ = score_numba_vectorized(
    np.ones(warmup_size, dtype=np.float64),  # reads_dna_normal
    np.ones(warmup_size, dtype=np.float64),  # reads_dna_tumor
    np.ones(warmup_size, dtype=np.float64),  # reads_rna_normal
    np.ones(warmup_size, dtype=np.float64),  # reads_rna_tumor
    np.ones(warmup_size, dtype=np.float64),  # reads_onekg
    np.ones(warmup_size, dtype=np.float64),  # samples_dna_normal
    np.ones(warmup_size, dtype=np.float64),  # samples_dna_tumor
    np.ones(warmup_size, dtype=np.float64),  # samples_rna_normal
    np.ones(warmup_size, dtype=np.float64),  # samples_rna_tumor
    np.ones(warmup_size, dtype=np.float64),  # samples_onekg
    np.ones(warmup_size, dtype=np.float64),  # pop_size_dna_normal
    np.ones(warmup_size, dtype=np.float64),  # pop_size_dna_tumor
    np.ones(warmup_size, dtype=np.float64),  # pop_size_rna_normal
    np.ones(warmup_size, dtype=np.float64),  # pop_size_rna_tumor
    np.full(warmup_size, 1.0),               # pop_size_dna_onekg
    np.full(warmup_size, 0.5),               # w_tumor
    np.full(warmup_size, 0.5),               # w_dna
    np.full(warmup_size, 0.5),               # w_read
    np.full(warmup_size, 50.0)               # upper_factor
)
print(f"Numba warmup complete (using {args.n_threads} threads)")
    
batch_size = args.batch_size # process in batches to save memory
for i in range(0, total_rows, batch_size):
    t_0 = time.time()
    batch_num = i // batch_size + 1
    total_batches = (total_rows + batch_size - 1) // batch_size
    print(f"Processing batch {batch_num}/{total_batches} (rows {i} to {min(i + batch_size, total_rows)})")

    df = conn.execute(f'SELECT * FROM read_csv_auto("{args.input}", delim="\t") LIMIT {batch_size} OFFSET {i}').df()
    
    # Helper to convert scalars to arrays (numba requires all inputs to be same-length arrays)
    def get_col_or_scalar(param_val, df_obj):
        if isinstance(param_val, str):
            return df_obj[param_val].values
        else:
            # Broadcast scalar to array of same length as dataframe
            return np.full(len(df_obj), param_val, dtype=np.float64)
    
    df['fusion_score'] = score_numba_vectorized(
        # reads
        get_col_or_scalar(params['read_dna_normal'], df),
        get_col_or_scalar(params['read_dna_tumor'], df),
        get_col_or_scalar(params['read_rna_normal'], df),
        get_col_or_scalar(params['read_rna_tumor'], df),
        get_col_or_scalar(params['read_dna_onekg'], df),
        # samples
        get_col_or_scalar(params['sample_dna_normal'], df),
        get_col_or_scalar(params['sample_dna_tumor'], df),
        get_col_or_scalar(params['sample_rna_normal'], df),
        get_col_or_scalar(params['sample_rna_tumor'], df),
        get_col_or_scalar(params['sample_dna_onekg'], df),
        # population sizes (broadcast scalars to arrays)
        get_col_or_scalar(params['pop_size_dna_normal'], df),
        get_col_or_scalar(params['pop_size_dna_tumor'], df),
        get_col_or_scalar(params['pop_size_rna_normal'], df),
        get_col_or_scalar(params['pop_size_rna_tumor'], df),
        get_col_or_scalar(params['pop_size_dna_onekg'], df),
        # hyperparameters (broadcast scalars to arrays)
        get_col_or_scalar(params['weight_tumor'], df),
        get_col_or_scalar(params['weight_dna'], df),
        get_col_or_scalar(params['weight_read'], df),
        get_col_or_scalar(params['upper_factor'], df)
    )
        
    
    # Write batch to TSV
    print(f"Writing batch {batch_num} to {args.output}")
    mode = "w" if i == 0 else "a"
    header = i == 0
    df.to_csv(args.output, sep="\t", index=False, mode=mode, header=header)
    print(f"Batch {batch_num} completed")
    print(f"Time for batch {batch_num}: {time.time() - t_0:.2f} seconds")

print(f"Wrote all batches to {args.output}")