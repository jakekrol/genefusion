#!/usr/bin/env python

import pandas as pd
import argparse
import os
parser=argparse.ArgumentParser()
parser.add_argument('-k', type=int, default=3)
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)

def filename2svtype(x):
    y= x.split('.')[0].split('_')[-1]
    if y == 'inversions':
        return 'INV'
    elif y == 'del':
         return 'DEL'
    elif y == 'dup':
         return 'DUP'
    elif y == 'translocations':
        return 'TRA'
    else:
        return 'BND'

       

args=parser.parse_args()

df = pd.read_csv(args.input, sep='\t', skiprows=1)
df['svtype'] = filename2svtype(os.path.basename(args.input))
df['tumor_paired_plus_split_samplewise_evidence'] = df['Pairend'] + df['Split']
df=df.sort_values('tumor_paired_plus_split_samplewise_evidence', ascending=False).reset_index(drop=True)
df_topk = df.head(args.k)
df_topk.to_csv(args.output, sep='\t', index=False)
