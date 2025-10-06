#!/usr/bin/env python3
import duckdb
import argparse
import numpy as np
from genefusion.genefusion import *
import yaml
import os,sys

parser = argparse.ArgumentParser(description='Score fusions')
parser.add_argument('-i', '--input', type=str,required=True, help='Input fusion feature file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output scored fusion file')
parser.add_argument('--w_dna', type=float, default=0.5, help='Weight for DNA evidence in scoring (default: 0.5)')
parser.add_argument('--upper_factor', type=float, default=100, help='Upper factor for scoring (default: 100)')
parser.add_argument('--pop_size_yaml',type=str,help='YAML file with population sizes (default: None)', default=None)
parser.add_argument('--col_map_yaml',type=str,help='Column mapping for input file',
                    default='../../config/score_fusion_colmap.yaml')
args = parser.parse_args()

conn = duckdb.connect()
total_rows = conn.execute(f'SELECT COUNT(*) FROM read_csv_auto("{args.input}", delim="\t")').fetchone()[0]
print(f"Total rows: {total_rows}")


with open(args.col_map_yaml) as f:
    colmap = yaml.safe_load(f)

batch_size = 500000  # Process in smaller batches
for i in range(0, total_rows, batch_size):
    batch_num = i // batch_size + 1
    total_batches = (total_rows + batch_size - 1) // batch_size
    print(f"Processing batch {batch_num}/{total_batches} (rows {i} to {min(i + batch_size, total_rows)})")

    df = conn.execute(f'SELECT * FROM read_csv_auto("{args.input}", delim="\t") LIMIT {batch_size} OFFSET {i}').df()
    df["fusion_score"] = np.vectorize(score)(
        # reads
        (df[colmap['read_dna_normal']] if colmap['read_dna_normal'] != str(0) else 0),
        (df[colmap['read_dna_tumor']] if colmap['read_dna_tumor'] != str(0) else 0),
        (df[colmap['read_rna_normal']] if colmap['read_rna_normal'] != str(0) else 0),
        (df[colmap['read_rna_tumor']] if colmap['read_rna_tumor'] != str(0) else 0),
        (df[colmap['read_onekg_dna']] if colmap['read_onekg_dna'] != str(0) else 0),
        # samples
        (df[colmap['sample_dna_normal']] if colmap['sample_dna_normal'] != str(0) else 0),
        (df[colmap['sample_dna_tumor']] if colmap['sample_dna_tumor'] != str(0) else 0),
        (df[colmap['sample_rna_normal']] if colmap['sample_rna_normal'] != str(0) else 0),
        (df[colmap['sample_rna_tumor']] if colmap['sample_rna_tumor'] != str(0) else 0),
        (df[colmap['sample_onekg_dna']] if colmap['sample_onekg_dna'] != str(0) else 0),
        # population sizes
        67, 70, 0, 70, 2536,
        # hyperparameters
        w_dna=args.w_dna,
        upper_factor=args.upper_factor,
    )
    
    # write for first batch, append for others
    mode = "w" if i == 0 else "a"
    header = i == 0
    df.to_csv(args.output, sep="\t", index=False, mode=mode, header=header)
    print(f"Batch {batch_num} completed")

print(f"Wrote all batches to {args.output}")