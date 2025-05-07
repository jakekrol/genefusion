#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import argparse
import re

# input: icgc metadata file
# output: file id to sample type mapping table

parser = argparse.ArgumentParser(description='Build a file id to sample type mapping')
parser.add_argument('-m', '--metadata_file', type=str, default='/data/jake/genefusion/data/meta/icgc25k-legacy-data-locations.tsv')
parser.add_argument('-o', '--output', type=str, default='/data/jake/genefusion/data/meta/fileid2sampletype.tsv')
args = parser.parse_args()

df = pd.read_csv(args.metadata_file, sep='\t', index_col=0, low_memory=False)
# filter out non-BAM files
df_bam = df[df['Format'] == 'BAM'].copy()
df_bam['Specimen_Type'] = df_bam['Specimen_Type'].str.lower()
df_bam.reset_index(drop=True, inplace=True)
X = []
n = df_bam.shape[0]
for i,row in df_bam.iterrows():
    file_id = row['File_ID']
    specimen = row['Specimen_Type']
    if pd.isna(specimen):
        print(f'File ID: {file_id} Specimen Type: {specimen} not recognised')
        continue
    if 'tumour' in specimen:
        specimen = 'tumour'
    elif 'normal' in specimen:
        specimen = 'normal'
    else:
        raise ValueError(f'Specimen type {specimen} not recognised')
    X.append([file_id, specimen])
    print(f'Processed {i+1}/{n} rows')
X = pd.DataFrame(X, columns=['File_ID', 'Specimen_Type'])
X.to_csv(args.output, sep='\t', index=False)