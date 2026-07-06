#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from polymerization.score import *
import yaml
import multiprocessing as mp

parser = argparse.ArgumentParser(description='Calculate fusion scores')
parser.add_argument('--inputdir', '-i', default='chunks')
parser.add_argument('--outdir', '-o', default='chunks_scored')
parser.add_argument('--colmap', default='score_column_map.yaml')
parser.add_argument('--cpus', '-c', type=int, default=30, help='number of cpus to use for parallel processing')
parser.add_argument('--chunk_columns', default = 'chunk_columns.txt')

def score_chunk(chunk_file, outfile, tissue, column_map_tissue, chunk_columns):
	# read
	df = pd.read_csv(chunk_file, sep='\t')
	df.columns = chunk_columns
	cols_keep = [
		'gene_left',
		'gene_right',
		'reads_high_coverage_1000g_dna',
		'samples_high_coverage_1000g_dna'
	]
	# subset cols based on tissue
	for col in df.columns.tolist():
		if tissue in col:
			cols_keep.append(col)
	df = df[cols_keep]
	# apply normalization
	df_norm = normalize_evidence_columns(df, column_map_tissue)
	# scale read and samples columns to [0,0.5]
	for col in df_norm.columns.tolist():
		if (col != 'gene_left') and col != ('gene_right'):
			df_norm[col] = df_norm[col] * 0.5
	# count num tumor and normal modalities
	n_tumor=0
	n_normal=0
	normal_cols = []
	tumor_cols = []
	for k,sub_dict in column_map_tissue.items():
		evidence_type = sub_dict['evidence_type']
		specimen = sub_dict['specimen']
		# only increment over sample (not read too) evidence to avoid double counting
		if evidence_type == 'sample':
			if specimen == 'tumor':
				n_tumor += 1
			elif specimen == 'normal':
				n_normal += 1
		# track tumor and normal columns for scoring
		if specimen == 'tumor':
			tumor_cols.append(k)
		elif specimen == 'normal':
			normal_cols.append(k)
		else:
			print(f"Warning: specimen type {specimen} not recognized for column {k}")
	# apply weight and negatively weight normal columns
	for k, sub_dict in column_map_tissue.items():
		specimen = sub_dict['specimen']
		if specimen == 'tumor':
			df_norm[k] = df_norm[k] * (1 / n_tumor)
		elif specimen == 'normal':
			df_norm[k] = df_norm[k] * (1 / n_normal) * -1
	# compute final score
	df_norm['score_uniform'] = df_norm[column_map_tissue.keys()].sum(axis=1)
	df_norm = df_norm.sort_values(by='score_uniform', ascending=False).reset_index(drop=True)
	# subset cols for output
	cols_out = ['gene_left', 'gene_right', 'score_uniform']
	df_norm = df_norm[cols_out]
	df_norm.to_csv(outfile, sep='\t', index=False)


def main():
	TISSUES=[
		'blood',
		'bone',
		'breast',
		'esophagus',
		'gallbladder',
		'headneck',
		'kidney',
		'liver',
		'ovary',
		'pancreas',
		'prostate'
	]
	args = parser.parse_args()
	os.makedirs(args.outdir, exist_ok=True)

	# inputs
	with open(args.colmap, 'r') as f:
		column_map = yaml.safe_load(f)
	chunkfiles = [os.path.join(args.inputdir, f) for f in os.listdir(args.inputdir)]
	chunkfiles.sort()
	chunk_columns = []
	with open(args.chunk_columns, 'r') as f:
		for line in f:
			chunk_columns.append(line.strip())
	for tissue in TISSUES:
		os.makedirs(os.path.join(args.outdir, tissue), exist_ok=True)
		column_map_tissue = {k: v for k, v in column_map.items() if (tissue in k) or (k == 'reads_high_coverage_1000g_dna') or (k == 'samples_high_coverage_1000g_dna')}
		arguments = []
		for chunk_file in chunkfiles:
			outfile = os.path.join(args.outdir, tissue, f"{os.path.basename(chunk_file)}_score.tsv")
			arguments.append((chunk_file, outfile, tissue, column_map_tissue, chunk_columns))
		with mp.Pool(processes=args.cpus) as pool:
			pool.starmap(score_chunk, arguments)

if __name__ == '__main__':
    main()