#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
from polymerization.giggle2fusion import *
from polymerization.stix2fusion import *
from polymerization.io import *
from polymerization.polymerization import *
import time

parser = argparse.ArgumentParser(description='Get gene fusion evidence 1000 Genomes high coverage data')
parser.add_argument('--shardfile', default='./shardfile.giggle.tsv', help='path to giggle shardfile')
parser.add_argument('--bedfile', default='../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed')
parser.add_argument('--fusion_set', default='../2025_11-fusion_lexicographic_key/gene_pairs_sorted_by_genomic_position.parquet')
parser.add_argument('--outdir', default='./g2f_out', help='output directory')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

### read inputs
print(f"# reading giggle shard file: {args.shardfile}")
df_shard = read_giggle_shardfile(args.shardfile)
print(f"# reading bed file: {args.bedfile}")
df_bed = read_bed(args.bedfile, gene_col_idx=3)

print(f"# reading fusion set: {args.fusion_set}")
t_0=time.time()
df_fusions = read_fusion_set(args.fusion_set, reader='polars')
print(f"# time to read fusion set: {time.time() - t_0:.2f} seconds")

# already sorted
# print(f"# sorting fusion set by genomic position")
# df_fusions_sort = left_sort_fusion_set(df_fusions, df_bed)

print(f"# converting fusion set and bed file to polars dataframes")
t_0=time.time()
df_fusions_pl = pl.from_pandas(df_fusions)
df_bed_pl = pl.from_pandas(df_bed)
print(f"# time to convert fusion set and bed file to polars dataframes: {time.time() - t_0:.2f} seconds")

### merge fusion set with bed file
print(f"# merging fusion set with bed file")
# set columns in fusion set
df_fusions_pl.columns = ['gene_left', 'gene_right']
t_0=time.time()
df_merge = merge_fusion_set_with_bed(df_fusions_pl, df_bed_pl, merger='polars')
print(f"# time to merge fusion set with bed file: {time.time() - t_0:.2f} seconds")

print(f"# converting merged dataframe back to pandas dataframe")
t_0=time.time()
df_merge_pd = df_merge.to_pandas()
print(f"# time to convert merged dataframe back to pandas dataframe: {time.time() - t_0:.2f} seconds")
print(f"# shape of merged dataframe: {df_merge_pd.shape}")

df_giggle = merge_fusion_set_bed2giggle(
    df_merge_pd, df_shard, outdir=args.outdir, max_workers=30, timeout=60*60*2, bgzip=True
)

df_giggle.to_csv(os.path.join(args.outdir, 'g2f_giggle.tsv'), sep='\t', index=False)
