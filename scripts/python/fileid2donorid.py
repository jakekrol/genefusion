#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import argparse
import re

# input: icgc file id
# output: donor id

### meta
df_meta = pd.read_csv('/data/jake/genefusion/data/meta/icgc25k-legacy-data-locations.tsv', sep='\t',index_col=0, low_memory=False)  

parser = argparse.ArgumentParser(description='Map file id to sample type {tumour/normal}')
parser.add_argument('-f', '--file_id', type=str, required=True)
args = parser.parse_args()

# search
x = df_meta[df_meta['File_ID'] == args.file_id]#['Specimen_Type'].str.lower().values
assert x.shape[0] == 1, f'File ID " {args.file_id} " search result is not singular. Shape is {x.shape}'
# get specimen type (tumour/normal)
x = x['ICGC_Donor'].values[0]
print(x)