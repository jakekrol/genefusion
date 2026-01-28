#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time
import subprocess
parser = argparse.ArgumentParser(description='Compute the total burden of each gene')
# burden is sum of all read supp. for any fusion involving the gene
parser.add_argument('-i', '--input', type=str, required=True, help='Input paired-end count file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for gene burden table')
parser.add_argument('--header', action='store_true', help='Input file has header')
parser.add_argument('--burden_col_name', type=str, default='burden_total', help='Name of the burden column in output file')
args = parser.parse_args()
burden = {}
with open(args.input, 'r') as f:
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
# sort burden by counts descending
burden = dict(sorted(burden.items(), key=lambda item: item[1], reverse=True))
with open(args.output, 'w') as out:
    out.write(f'gene\t{args.burden_col_name}\n')
    for gene, burden in burden.items():
        out.write(f'{gene}\t{burden}\n')
