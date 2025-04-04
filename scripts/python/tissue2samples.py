#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import argparse
import re

### meta
df_abbrevs = pd.read_csv('/data/jake/genefusion/data/meta/icgc-abbrevs-simple.csv')
df_meta = pd.read_csv('/data/jake/genefusion/data/meta/icgc25k-legacy-data-locations.tsv', sep='\t',index_col=0, low_memory=False)  
df_meta = df_meta[df_meta['Study']== 'PCAWG']
df_meta['Project'] = df_meta['Project'].str.split('-').str[0]
### args
parser = argparse.ArgumentParser(description='Get sample counts for PCAWG tissues')
parser.add_argument('-t', '--tissue', type=str, required=True)
parser.add_argument('-f', '--outformat', type=str, required=False, default='table')
parser.add_argument('-o', '--output', type=str, required=False, default='pcawg_meta.tsv')
parser.add_argument('-c', '--cancer_only', required=False, action='store_true')
parser.add_argument('-n', '--normal_only', required=False, action='store_true')
parser.add_argument('-s', '--seqtype', type=str, required=True)

args = parser.parse_args()

formats = set(['table', 'file_ids', 'counts'])
if args.outformat not in formats:
    raise ValueError(f'Output format must be one of {formats}')

if (args.cancer_only or args.normal_only) and args.outformat == 'counts':
    raise ValueError('Counts cannot be used with cancer or normal only')

seqtypes = set(['DNA', 'RNA'])
if args.seqtype not in seqtypes:
    raise ValueError(f'Sequencing type must be one of {seqtypes}')
    

# get projects
projects = set(df_abbrevs[df_abbrevs['tissue'] == args.tissue]['abbreviation'])

# subset meta
df_meta = df_meta[df_meta['Project'].isin(projects)]

# only consider bams
df_meta = df_meta[df_meta['Data_Type'] == 'BAM']
# drop minibams
df_meta = df_meta[~df_meta['File_Name'].str.contains('mini')]

# seqtype
if args.seqtype == 'DNA':
    df_meta = df_meta[df_meta['Experimental_Strategy'].str.contains('WGS')]
if args.seqtype == 'RNA':
    df_meta = df_meta[df_meta['Experimental_Strategy'].str.contains('RNA-Seq')]

if args.cancer_only:
    df_meta['Specimen_Type'] = df_meta['Specimen_Type'].str.lower()
    df_meta = df_meta[df_meta['Specimen_Type'].str.contains('tumour')]
        
if args.normal_only:
    df_meta['Specimen_Type'] = df_meta['Specimen_Type'].str.lower()
    df_meta = df_meta[df_meta['Specimen_Type'].str.contains('normal')]

### table
if args.outformat == 'table':
    df_meta.to_csv(args.output, sep='\t', index=False)
    sys.exit()

### counts
if args.outformat == 'counts':
    df_meta['Specimen_Type'] = df_meta['Specimen_Type'].str.lower()
    n_cancer = df_meta[df_meta['Specimen_Type'].str.contains('tumour')].shape[0]
    n_normal = df_meta[df_meta['Specimen_Type'].str.contains('normal')].shape[0]
    print(f'Cancer: {n_cancer}')
    print(f'Normal: {n_normal}')
    print(f'Total: {n_cancer + n_normal}')
    sys.exit()

### file ids
file_ids = df_meta['File_ID'].tolist()
with open(args.output, 'w') as f:
    for file_id in file_ids:
        f.write(f'{file_id}\n')
    