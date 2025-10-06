#!/usr/bin/env bash
icgc_abbrev='icgc-cancer-abbreviations.txt'
script="import pandas as pd;\
icgc_abbrev='icgc-cancer-abbreviations.txt';\
df = pd.read_csv(icgc_abbrev, sep=':', skiprows=[0,1],header=None);\
df['site']=df[1].apply(lambda x: x.split('(')[1]);\
df['site']=df['site'].apply(lambda x: x.split(',')[0]);\
df.columns=['code','name','site'];\
df.to_csv('icgc_code_name_site.tsv', sep='\t', index=False);"
python -c "$script"





