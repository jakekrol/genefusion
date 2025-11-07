#!/usr/bin/env python3
import duckdb
import swifter
import dask
import argparse
import numpy as np
from genefusion.genefusion import score
import yaml
import os,sys
# from pandarallel import pandarallel

swifter.set_defaults(npartitions=os.cpu_count()-10 if os.cpu_count()>10 else 1)
swifter.set_defaults(progress_bar=True)
swifter.set_defaults(allow_dask_on_strings=True)





parser = argparse.ArgumentParser(description='Score fusions')
parser.add_argument('-i', '--input', type=str,required=True, help='Input fusion feature file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output scored fusion file')
parser.add_argument('--score_yaml',type=str,help='YAML file with population sizes and column names', required=True)
# parser.add_argument('--batch_size',type=int,default=20000000,help='Number of rows to process per batch (default: 10000000)')
parser.add_argument('--batch_size',type=int,default=20000000,help='Number of rows to process per batch (default: 20000000)')
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

batch_size = args.batch_size # process in batches to save memory
for i in range(0, total_rows, batch_size):
    batch_num = i // batch_size + 1
    total_batches = (total_rows + batch_size - 1) // batch_size
    print(f"Processing batch {batch_num}/{total_batches} (rows {i} to {min(i + batch_size, total_rows)})")

    df = conn.execute(f'SELECT * FROM read_csv_auto("{args.input}", delim="\t") LIMIT {batch_size} OFFSET {i}').df()
    df['fusion_score'] = df.swifter.apply(
    # df['fusion_score'] = df.parallel_apply(
        lambda row: score(
            # reads
            row[params['read_dna_normal']] if type(params['read_dna_normal']) == str else int(params['read_dna_normal']),
            row[params['read_dna_tumor']] if type(params['read_dna_tumor']) == str else int(params['read_dna_tumor']),
            row[params['read_rna_normal']] if type(params['read_rna_normal']) == str else int(params['read_rna_normal']),
            row[params['read_rna_tumor']] if type(params['read_rna_tumor']) == str else int(params['read_rna_tumor']),
            row[params['read_dna_onekg']] if type(params['read_dna_onekg']) == str else int(params['read_dna_onekg']),
            # samples
            row[params['sample_dna_normal']] if type(params['sample_dna_normal']) == str else int(params['sample_dna_normal']),
            row[params['sample_dna_tumor']] if type(params['sample_dna_tumor']) == str else int(params['sample_dna_tumor']),
            row[params['sample_rna_normal']] if type(params['sample_rna_normal']) == str else int(params['sample_rna_normal']),
            row[params['sample_rna_tumor']] if type(params['sample_rna_tumor']) == str else int(params['sample_rna_tumor']),
            row[params['sample_dna_onekg']] if type(params['sample_dna_onekg']) == str else int(params['sample_dna_onekg']),
            # population sizes
            params['pop_size_dna_normal'],
            params['pop_size_dna_tumor'],
            params['pop_size_rna_normal'],
            params['pop_size_rna_tumor'],
            params['pop_size_dna_onekg'],
            # hyperparameters
            w_tumor=params['weight_tumor'],
            w_dna=params['weight_dna'],
            w_read=params['weight_read'],
            upper_factor=params['upper_factor'],
        ),
        axis=1
    )
    # df["fusion_score"] = np.vectorize(score)(
    #     # reads
    #     (df[params['read_dna_normal']] if type(params['read_dna_normal']) == str else int(params['read_dna_normal'])),
    #     (df[params['read_dna_tumor']] if type(params['read_dna_tumor']) == str else int(params['read_dna_tumor'])),
    #     (df[params['read_rna_normal']] if type(params['read_rna_normal']) == str else int(params['read_rna_normal'])),
    #     (df[params['read_rna_tumor']] if type(params['read_rna_tumor']) == str else int(params['read_rna_tumor'])),
    #     (df[params['read_dna_onekg']] if type(params['read_dna_onekg']) == str else int(params['read_dna_onekg'])),
    #     # samples
    #     (df[params['sample_dna_normal']] if type(params['sample_dna_normal']) == str else int(params['sample_dna_normal'])),
    #     (df[params['sample_dna_tumor']] if type(params['sample_dna_tumor']) == str else int(params['sample_dna_tumor'])),
    #     (df[params['sample_rna_normal']] if type(params['sample_rna_normal']) == str else int(params['sample_rna_normal'])),
    #     (df[params['sample_rna_tumor']] if type(params['sample_rna_tumor']) == str else int(params['sample_rna_tumor'])),
    #     (df[params['sample_dna_onekg']] if type(params['sample_dna_onekg']) == str else int(params['sample_dna_onekg'])),
    #     # population sizes
    #     params['pop_size_dna_normal'],
    #     params['pop_size_dna_tumor'],
    #     params['pop_size_rna_normal'],
    #     params['pop_size_rna_tumor'],
    #     params['pop_size_dna_onekg'],
    #     # hyperparameters
    #     w_tumor=params['weight_tumor'],
    #     w_dna=params['weight_dna'],
    #     w_read=params['weight_read'],
    #     upper_factor=params['upper_factor'],
    # )
    
    # write for first batch, append for others
    mode = "w" if i == 0 else "a"
    header = i == 0
    print(f"Writing batch {batch_num} to {args.output}")
    df.to_csv(args.output, sep="\t", index=False, mode=mode, header=header)
    print(f"Batch {batch_num} completed")

print(f"Wrote all batches to {args.output}")