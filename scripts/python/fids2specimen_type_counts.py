#!/usr/bin/env python3
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="Count specimen types from file IDs")
parser.add_argument("--fileids", required=True, help="File containing line-delim list of pcawg file IDs")
parser.add_argument("--lookup", required=False, help="TSV file mapping File_ID to Specimen_Type",
                    default='/data/jake/genefusion/results/2025_05-pcawg_fileid2sample_type/fileid2sampletype.tsv')
args = parser.parse_args()
df_l=pd.read_csv(args.lookup,sep='\t')
with open(args.fileids, 'r') as file:
    x = file.read().splitlines()
mask = df_l.File_ID.isin(x)
df_lm=df_l[mask]
print(df_lm['Specimen_Type'].value_counts())
