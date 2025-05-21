#!/usr/bin/env python3
import os
import sys
import argparse
import shutil
import multiprocessing as mp

parser = argparse.ArgumentParser(description='Migrate files by specimen type')
parser.add_argument('-i', '--input', type=str, required=True, help='Input directory with files to migrate')
parser.add_argument('-t', '--dir_tumour', type=str, required=True, help='Tumour directory to migrate files to')
parser.add_argument('-n', '--dir_normal', type=str, required=True, help='Normal directory to migrate files to')
parser.add_argument('-p', '--processes', type=int, default=mp.cpu_count(), help='Number of processes to use for migration')
args = parser.parse_args()

def fname2type(fname):
    if 'tumour' in fname.lower():
        return 'tumour'
    elif 'normal' in fname.lower():
        return 'normal'
    else:
        raise ValueError(f"Unknown specimen type for file: {fname}")
files = os.listdir(args.input)
files = [os.path.join(args.input, f) for f in files]
locations = []
for src in files:
    specimen_type = fname2type(src)
    # remove the prefix from the filename
    dst = os.path.basename(src)
    # remove the specimen type from the filename
    dst = dst.replace('tumour.', '').replace('normal.', '')
    if specimen_type == 'tumour':
        dst = os.path.join(args.dir_tumour, dst)
    elif specimen_type == 'normal':
        dst = os.path.join(args.dir_normal, dst)
    else:
        raise ValueError(f"Unknown specimen type for file: {src}")
    locations.append((src, dst))

print(f'Migrating {len(locations)} files to {args.dir_tumour} and {args.dir_normal} using {args.processes} processes')
with mp.Pool(args.processes) as pool:
    pool.starmap(shutil.copy, locations)


    
