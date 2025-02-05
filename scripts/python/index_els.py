#!/usr/bin/env python3

# sub names in edge lists with index

import pandas as pd
import os,sys
import time

INDEX = '/data/jake/genefusion/data/genes.index'

EXAMPLE='NA12347.fusion'

df_idx = pd.read_csv(INDEX, sep='\t', header=None, names=['gene'], usecols = [1])
print(df_idx.head())

df_el = pd.read_csv(EXAMPLE, sep='\t', header=None, names=['gene1', 'gene2'])
print(df_el.head())

z = 2
t_s=time.time()
with open('test2.txt', 'w') as f:
    for i,row in df_el.iterrows():
        gene1 = row['gene1'].split('.')[0] # remove version number
        gene2 = row['gene2'].split('.')[0] # remove version number
        idx1 = str(df_idx[df_idx['gene'] == gene1].index[0])
        idx2 = str(df_idx[df_idx['gene'] == gene2].index[0])
        f.write(idx1 + '\t' + idx2 + '\n')
# for i,row in df_el.iterrows():
    # gene1 = row['gene1'].split('.')[0] # remove version number
    # gene2 = row['gene2'].split('.')[0] # remove version number
    # idx1 = df_idx[df_idx['gene'] == gene1].index
    # idx2 = df_idx[df_idx['gene'] == gene2].index[0]
    # df_el_idx = pd.concat([df_el_idx, pd.DataFrame({'gene1': idx1, 'gene2': idx2})], ignore_index=True)
t_e = time.time()
print('time:', t_e-t_s)
# df_el_idx = df_el_idx.astype(int)
