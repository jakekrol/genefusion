#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time
import subprocess

parser = argparse.ArgumentParser(description='Subset fusion table for query')
parser.add_argument('-i', '--input', type=str, required=True, help='Input gene fusion table file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file')  
parser.add_argument('-s', '--script_path', type=str, default='/data/jake/genefusion/scripts/python/leftgene.py', help='Path to the leftgene script')
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=None)
assert df.shape[1] == 2, f'Input file {args.input} should have 2 columns'

# get the leftgenes
def leftgene(x,y):
    result = subprocess.run(['python3', args.script_path, '-x', x, '-y', y], capture_output=True, text=True)
    out = result.stdout.strip()
    return out
with open(args.output, 'w') as out:
    out.write('left\tright\tsort_result\n')
    n = df.shape[0]
    i = 1
    for _,row in df.iterrows():
        x = row[0]
        y = row[1]
        result = leftgene(x,y)
        if result == "-1":
            print(f'{x} or {y} not found in the gene bed file')
            out.write(f'{x}\t{y}\t-1\n')
        else:
            print(f'leftgene: {x} and {y} -> {result}')
            if x == result:
                out.write(f'{x}\t{y}\t1\n')
            elif y == result:
                out.write(f'{y}\t{x}\t1\n')
            else:
                print(f'unexpected leftgene result: "{result}" for "{x}" and "{y}"')
                out.write(f'{x}\t{y}\t-1\n')
        print(f'Processed {i}/{n} rows')
        i+=1
