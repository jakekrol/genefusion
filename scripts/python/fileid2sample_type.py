#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import argparse
import re

# input: icgc file id
# output: tumour/normal

### meta
# df_abbrevs = pd.read_csv('/data/jake/genefusion/data/meta/icgc-abbrevs-simple.csv')
df_meta = pd.read_csv('/data/jake/genefusion/data/meta/icgc25k-legacy-data-locations.tsv', sep='\t',index_col=0, low_memory=False)  
# df_meta = df_meta[df_meta['Study']== 'PCAWG']
# df_meta['Project'] = df_meta['Project'].str.split('-').str[0]

parser = argparse.ArgumentParser(description='Map file id to sample type {tumour/normal}')
parser.add_argument('-f', '--file_id', type=str, required=True)
args = parser.parse_args()

# search
x = df_meta[df_meta['File_ID'] == args.file_id]#['Specimen_Type'].str.lower().values
assert x.shape[0] == 1, f'File ID " {args.file_id} " search result is not singular. Shape is {x.shape}'
# check format is BAM
fmat = x['Format'].values[0]
assert fmat == 'BAM', f'File ID " {args.file_id} " not a BAM file'
# get specimen type (tumour/normal)
x = x['Specimen_Type'].values[0]
if 'tumour' in x.lower():
    print('tumour')
elif 'normal' in x.lower():
    print('normal')
else:
    raise ValueError(f'Specimen type {x} not recognised')