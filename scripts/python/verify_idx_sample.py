#!/usr/bin/env python3
import os,sys
import glob
import pandas as pd
import argparse
from collections import Counter

# input1: dir with giggle index file ids
# input2: file with the corresponding tissue's DNA tumour samples
# input3: file with the corresponding tissue's DNA normal samples

parser = argparse.ArgumentParser(description='Verify the samples in the giggle index are tumour samples')
parser.add_argument('--giggle_dir', '-g', type=str, help='Directory with giggle index file ids')
parser.add_argument('--tumour_samples', '-t', type=str, help='File with the corresponding tissue\'s DNA tumour samples')
parser.add_argument('--normal_samples', '-n', type=str, help='File with the corresponding tissue\'s DNA normal samples')
parser.add_argument('--output', '-o', type=str, help='Report number/fraction of index samples that are tumour/normal and which')
args = parser.parse_args()

# giggle index samples
G = os.listdir(args.giggle_dir)
G = [os.path.basename(x).split('.')[0] for x in G]
# check for dups
C = dict(Counter(G))
for v in C.values():
    if v > 1:
        print(f'Warning: duplicate sample ID found in index: " {v} "')
print(C)
G = set(G)

# tumour samples
T = set()
with open(args.tumour_samples, 'r') as f:
    for line in f:
        sample = line.strip()
        T.add(sample)
# normal samples
N = set()
with open(args.normal_samples, 'r') as f:
    for line in f:
        sample = line.strip()
        N.add(sample)
print(G,T,N)

# count overlap
Tint = T.intersection(G)
Nint = N.intersection(G)
# intersection over union
Tunion = T.union(G)
Nunion = N.union(G)
# instead of printing, write to file
with open(args.output, 'w') as f:
    f.write(f'Number samples in index: {len(G)}\n')
    f.write(f'Number samples in tumour file: {len(T)}\n')
    f.write(f'Number samples in normal file: {len(N)}\n')
    f.write('---\n')
    f.write(f'Num index samples in tumour file: {len(Tint)}\n')
    f.write(f'Num index samples in normal file: {len(Nint)}\n')
    f.write(f'Fraction index samples in tumour file: {len(Tint) / len(G)}\n')
    f.write(f'Fraction index samples in normal file: {len(Nint) / len(G)}\n')
    # write which samples in index are tumour and which are normal
    f.write('---\n')
    f.write('Samples in index that are tumour:\n')
    for sample in Tint:
        f.write(f'{sample}\n')
    f.write('---\n')
    f.write('Samples in index that are normal:\n')
    for sample in Nint:
        f.write(f'{sample}\n')
    f.write('---\n')
    f.write('Samples in index that are neither:\n')
    for sample in G:
        if (sample not in Tint) and (sample not in Nint):
            f.write(f'{sample}\n')
    f.write('---\n')

sys.exit()


