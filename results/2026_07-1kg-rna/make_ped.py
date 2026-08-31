#!/usr/bin/env python

import argparse
import glob
from jkbiolib.datasets.loaders import *
import os
import pandas
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("--metadata", default="../../data/2026_08-thousg-rna-metadata/thousg-mage-sample-metadata.tsv")
parser.add_argument("--dir_bed", default="thousg_rna_bed_sort_cln")
parser.add_argument("--output", default="thousg_rna.ped")
parser.add_argument("--rename", action='store_true')
args = parser.parse_args()

def ftp_url2srr_id(url):
	return os.path.basename(os.path.dirname(url))

def rename_bed(path, srr2sample):
	directory = os.path.dirname(path)
	path = os.path.basename(path)
	srrid = path.split('.')[0]
	sample = srr2sample.get(srrid)
	path = f'{sample}.{path}'
	return os.path.join(directory, path)

def main():
	assert os.path.exists(args.metadata)
	assert os.path.isdir(args.dir_bed)
	# input
	df_meta = pd.read_csv(args.metadata, sep="\t")
	df_samples = thousg_rna_short_read_samples()
	df_samples['srr_id'] = df_samples['url'].apply(ftp_url2srr_id)
	# clean sample columns
	df_meta['Sample name'] = df_meta['Sample name'].apply(lambda x: x.strip())
	df_samples['Sample'] = df_samples['Sample'].apply(lambda x: x.strip())
	bed_files = glob.glob(os.path.join(args.dir_bed, "*.bed.gz"))
	srr2sample = dict(zip(df_samples["srr_id"], df_samples["Sample"]))
	# sub srrid to sample id in bed file name
	if args.rename:
		for f in bed_files:
			shutil.move(f, rename_bed(f, srr2sample))
	alt_files = glob.glob(os.path.join(args.dir_bed, "*.bed.gz"))
	# select relevant metadata
	cols_samples = ['Sample', 'srr_id', 'Data collection', 'Population']
	cols_meta =['Sample name', 'Sex', 'Biosample ID', 'Population code', 'Population name', 'Superpopulation code', 'Superpopulation name']
	df_samples = df_samples[cols_samples]
	df_samples.columns = map(lambda x: x.lower().replace(' ', '_'), df_samples.columns.tolist())
	df_meta = df_meta[cols_meta]
	df_meta.columns = map(lambda x: x.lower().replace(' ', '_'), df_meta.columns.tolist())
	df_merge = pd.merge(df_meta, df_samples, left_on='sample_name', right_on='sample', how = 'outer')
	na_values_metadata = df_merge.isna().sum().sum()
	print(f"# na_values_metadata={na_values_metadata}")
	order = ['sample', 'srr_id', 'biosample_id', 'sex', 'population_code', 'population_name', 'population', 'superpopulation_code', 'superpopulation_name', 'data_collection']
	df_merge = df_merge[order]
	df_merge['alt_file'] = ''
	n = len(alt_files)
	for i,f in enumerate(alt_files):
		if (i+1) % 100 == 0:
			print(f"# setting alt_file {i+1}/{n}")
		f=os.path.basename(f)
		sample = f.split('.')[0]
		srrid = f.split('.')[1]
		for j, row in df_merge.iterrows():
			if ((sample == row['sample']) and (srrid == row['srr_id'])):
				df_merge.loc[j, 'alt_file'] = f
	cols = df_merge.columns.tolist()
	cols.remove('alt_file')
	cols = ['alt_file'] + cols
	df_merge = df_merge[cols]
	df_merge = df_merge.reset_index(drop=True)
	mask = df_merge['alt_file'] == ''
	df_merge = df_merge[~mask]
	df_merge.drop_duplicates(inplace=True)
	df_merge.to_csv(args.output, index=False, sep="\t")

if __name__ == "__main__":
    main()