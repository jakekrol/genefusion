#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--icgc_metadata", type=str, default="../../data/2024_08-icgc_legacy_locations/icgc25k-legacy-data-locations.no_index.tsv")
parser.add_argument("--ega_datasets", type=str, default="../../data/2026_06-ega-datasets/datasets.txt")
parser.add_argument("--output_txt", default='summary.txt', type=str)
parser.add_argument("--output_tsv", default='icgc-ega-bams_fastqs.tsv', type=str)
parser.add_argument("--output_file_ids", default='icgc-ega-bams_fastqs.ids.txt', type=str)

args = parser.parse_args()

df = pd.read_csv(args.icgc_metadata, sep="\t")
datasets = set()
with open(args.ega_datasets, "r") as f:
	for line in f:
		datasets.add(line.strip())
mask = df['ega_dataset_id'].isin(datasets)
df_filtered = df[mask]
formats = {'BAM', 'FASTQ'}
mask = df_filtered['Format'].isin(formats)
df_filtered = df_filtered[mask]
df_filtered.to_csv(args.output_tsv,sep='\t',index=False)
with open(args.output_txt, "w") as f:
	x = df_filtered.groupby('Format').size()
	f.write("# file formats\n")
	for format, count in x.items():
		f.write(f"{format}: {count}\n")
	x = df_filtered.groupby('Experimental_Strategy').size()
	f.write("# experimental strategies\n")
	for strategy, count in x.items():
		f.write(f"{strategy}: {count}\n")

df_filtered = df_filtered.sort_values(by=['ega_file_id'])
with open(args.output_file_ids, "w") as f:
	for _, row in df_filtered.iterrows():
		f.write(f"{row['ega_file_id']}\n")
