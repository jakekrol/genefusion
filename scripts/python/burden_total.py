#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time
parser = argparse.ArgumentParser(description='Compute the total burden of each gene')
parser.add_argument('-i', '--input_fusion_table', type=str, required=True, help='Input fusion count table')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for gene burden table')
args = parser.parse_args()
df = pd.read_csv(args.input_fusion_table, sep='\t', header=None)
i = 1
burden={}
n = df.shape[0]
for _,row in df.iterrows():
    left = row[0]
    right = row[1]
    pe_count = row[2]
    if left != right:
        # add keys
        if left not in burden:
            burden[left] = 0
        if right not in burden:
            burden[right] = 0
        # update burden
        burden[left] += pe_count
        burden[right] += pe_count
    print(f'Processed {i}/{n} rows')
    i+= 1
# write the burden to the output file
with open(args.output, 'w') as out:
    for gene, burden in burden.items():
        out.write(f'{gene}\t{burden}\n')
