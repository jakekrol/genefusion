#!/usr/bin/env python3

# input fusion file
# output fusion files for each sample
# column 11 is sample

import pandas as pd
import os,sys

infile='/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/fusions_25_01_05/10.neg.A1CF.52559169.52645435.fusion'
outdir='/data/jake/genefusion/scratch/2025-01_17_25-sample_wise_fusions'

# group by sample and write to file 
df = pd.read_csv(infile, sep='\t', header=None)

gene = os.path.basename(infile).split('.')[2:-3][0]
print("Gene: ", gene) 

for sample, group in df.groupby(10):
    # strip the folder index prefix
    # and strip extension suffices
    sample = os.path.basename(sample)
    sample = sample.split('.')[0]

    group['source'] = gene
    group.rename(columns={3: 'target'}, inplace=True)
    # only bother keeping source and target gene column
    group = group[['source', 'target']]
    group.to_csv(f'{outdir}/{sample}.fusion', sep='\t', header=None, index=False)

# pd.read_csv(infile, sep='\t', header=None).\
#     groupby(10).apply(
#         lambda x: x.to_csv(f'{outdir}/{x}', sep='\t', header=None, index=False)
#     )