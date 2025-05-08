#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import argparse
import re

# input: icgc metadata file
# output: file name to file id mapping table

parser = argparse.ArgumentParser(description='Build a file name to file id mapping')
parser.add_argument('-m', '--metadata_file', type=str, default='/data/jake/genefusion/data/meta/icgc25k-legacy-data-locations.tsv')
parser.add_argument('-o', '--output', type=str, default='/data/jake/genefusion/data/meta/filename2fileid.tsv')
args = parser.parse_args()

df = pd.read_csv(args.metadata_file, sep='\t', index_col=0, low_memory=False)
# filter out non-BAM files
df_bam = df[df['Format'] == 'BAM'].copy()
# filter out mini BAMs
mask = df_bam['File_Name'].str.contains('mini')
df_bam = df_bam[~mask].copy()
# only keep PCAWG files
mask = df_bam['File_Name'].str.contains('PCAWG')
df_bam = df_bam[mask].copy()
df_bam.reset_index(drop=True, inplace=True)
X = []
n = df_bam.shape[0]
for i,row in df_bam.iterrows():
    fname = row['File_Name']
    file_id = row['File_ID']
    X.append([fname, file_id])
    print(f'Processed {i+1}/{n} rows')
X = pd.DataFrame(X, columns=['File_Name', 'File_ID'])
X.to_csv(args.output, sep='\t', index=False)