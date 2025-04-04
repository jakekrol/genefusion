#!/usr/bin/env python3
import sys,os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from genefusion.genefusion import cln_giggle
import subprocess
import time

bedtools='/data/jake/bedtools.static.binary'
genefile='/data/jake/genefusion/data/gene_file_cln.txt'
plotter='/data/jake/rl-tools/plot/ccdf.py'

# args
argparser = argparse.ArgumentParser(description='Query a fusion and output ccdf plot')
argparser.add_argument('-x','--gene_x', type=str, help='gene x')
argparser.add_argument('-y','--gene_y', type=str, help='gene y')
argparser.add_argument('-t','--tissue', type=str, help='tissue')
argparser.add_argument('-k','--keep', action='store_true', help='keep intermediate files')
argparser.add_argument('-o','--output', type=str, help='output base name')
argparser.add_argument('-i','--index', type=str, help='index directory')
argparser.add_argument('-b','--bed', type=str, help='gene location bed file')
args = argparser.parse_args()
if args.index is None:
    args.index = f'{args.tissue}_sort_b'
# check
if not os.path.exists(args.index):
    print(f'error: "{args.index}" does not exist')
    sys.exit(1)
# if args.bed is None:
#     args.bed = '/data/jake/genefusion/data/Homo_sapiens.GRCh37.82.genes.bed'

### giggle search
# choose gene with lexicographically less coordinates
df = pd.read_csv(genefile, sep='\t', header=None)
df.columns = ['chrm','left', 'right','gene','strand']
# lookup coords
try:
    # care of invalid unicode characters
    i = df.index[df['gene']==str(args.gene_x)][0]
except IndexError as e:
    print(e)
    print(f'error looking up gene " {repr(args.gene_x)} "  in bed file')
    sys.exit(1)
chr1 = df.loc[i,'chrm']
left1 = df.loc[i,'left']
right1 = df.loc[i,'right']
s1 = f'{chr1}{left1}{right1}'
try:
    j = df.index[df['gene']==args.gene_y][0]
except IndexError as e:
    print(e)
    print(f'error looking up gene " {repr(args.gene_y)} "  in bed file')
    sys.exit(1)
chr2 = df.loc[j,'chrm']
left2 = df.loc[j,'left']
right2 = df.loc[j,'right']
s2 = f'{chr2}{left2}{right2}'
# compare str representations
if s1 < s2:
    query=args.gene_x
    region=f'{chr1}:{left1}-{right1}'
else:
    query=args.gene_y
    region=f'{chr2}:{left2}-{right2}'
# check
out_giggle = f'{query}.giggle'
out_giggle_cln = f'{query}.cln.giggle'
if os.path.exists(out_giggle):
    print(f'warning: "{out_giggle}" already exists')
# search
print(f"query: {query},{region}")
print("search")
t=time.time()
cmd = f'giggle search -i {args.index} -r {region} -v > {out_giggle}' 
result = subprocess.run(cmd,shell=True,check=True)
print(f'search time: {time.time()-t}')
print(f'search result: {result}')
# clean giggle
print('clean')
df = pd.read_csv(out_giggle, sep='\t', header=None)
df = cln_giggle(df)
df.to_csv(out_giggle_cln, sep='\t', header=None, index=False)

### extract hits (right-hand read)
out_hits = f'{query}.hits'
print('cut')
cmd = f'cut -f 5- {out_giggle_cln} > {out_hits}'
result = subprocess.run(cmd,shell=True,check=True)
print(f'cut result: {result}')

### intersect
print("intersect")
out_intersect = f'{query}.intersect'
cmd = f'{bedtools} intersect -a {genefile} -b {out_hits} > {out_intersect}'
print(cmd)
result = subprocess.run(cmd,shell=True,check=True)
print(f'intersect result: {result}')

### count fusions
out_fusions = f'{query}.fusions'
cmd = f"cut -f 4 {out_intersect} | sort | uniq -c | sort -n > {out_fusions}"
result = subprocess.run(cmd,shell=True,check=True)
print(f'cut result: {result}')

### plot
genes = {args.gene_x,args.gene_y}
target = genes - {query}
target = target.pop()
# get count of target partner
cmd = f"cut -f 4 {out_intersect} | grep -c {target}"
try:
    result = subprocess.run(cmd,shell=True,check=True,capture_output=True,text=True)
    vline = int(result.stdout.strip()) # use this as vline
except:
    result = -1
    vline = int(0)
# get count distribution of all partners
cmd = f"sed -i 's/[a-zA-Z].*//g' {out_fusions}"
result = subprocess.run(cmd,shell=True,check=True)
print(f'sed result: {result}')
print("plot")
cmd = f"cat {out_fusions} | {plotter} -o p.png --xlog --title '{query} fusions ({target} marked)' -x 'Read evidence count' --axvline {vline}"
result = subprocess.run(cmd,shell=True,check=True)
print(f'plot result: {result}')
