#!/usr/bin/env python

import pandas as pd
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-k', type=int, default=3)
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)


args=parser.parse_args()

df = pd.read_csv(args.input, sep='\t', skiprows=1)
df['tumor_paired_plus_split_samplewise_evidence'] = df['Pairend'] + df['Split']
df=df.sort_values('tumor_paired_plus_split_samplewise_evidence', ascending=False).reset_index(drop=True)
df_topk = df.head(args.k)
df_topk.to_csv(args.output, sep='\t', index=False)
