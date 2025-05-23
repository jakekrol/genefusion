#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time
parser = argparse.ArgumentParser(description='Compute the total burden of each gene')
parser.add_argument('-i', '--input_pe_count', type=str, required=True, help='Input paired-end count file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for gene burden table')
parser.add_argument('--header', action='store_true', help='Input file has header')
args = parser.parse_args()
burden = {}
with open(args.input_pe_count, 'r') as f:
    if args.header:
        next(f)  # skip the header
    for i, line in enumerate(f):
        left, right, pe_count = line.strip().split('\t')
        pe_count = int(pe_count)
        left = str(left)
        right = str(right)
        if left != right:
            # add keys
            if left not in burden:
                burden[left] = 0
            if right not in burden:
                burden[right] = 0
            # update burden
            burden[left] += pe_count
            burden[right] += pe_count
        print(f'Processed {i+1} rows')
with open(args.output, 'w') as out:
    for gene, burden in burden.items():
        out.write(f'{gene}\t{burden}\n')
# if args.header:
#     df = pd.read_csv(args.input_pe_count, sep='\t')
# else:
#     df = pd.read_csv(args.input_pe_count, sep='\t', header=None)
# i = 1
# burden={}
# n = df.shape[0]
# for _,row in df.iterrows():
#     left = row[0]
#     right = row[1]
#     pe_count = row[2]
#     if left != right:
#         # add keys
#         if left not in burden:
#             burden[left] = 0
#         if right not in burden:
#             burden[right] = 0
#         # update burden
#         burden[left] += pe_count
#         burden[right] += pe_count
#     print(f'Processed {i}/{n} rows')
#     i+= 1
# # write the burden to the output file
# with open(args.output, 'w') as out:
#     for gene, burden in burden.items():
#         out.write(f'{gene}\t{burden}\n')
