#!/usr/bin/env python3

import os,sys
import argparse
import subprocess
import shutil
import time

# this script will take be slow until the cache is built
# consider pre-computing a lookup table of file ids to sample types before hand to parallelize this

# input: 1) a sample unswapped intersect file, 2) normal dir, 3) tumour dir, 4) path to file id to sample type script
# output: None 
# effect: copy file to normal/tumour directories

parser = argparse.ArgumentParser(description='Copy sample intersect file to normal/tumour directories')
parser.add_argument('-i', '--input_dir', type=str, required=True, help='Directory containing intersect files')
parser.add_argument('-n', '--normal_dir', type=str, required=True)
parser.add_argument('-t', '--tumour_dir', type=str, required=True)
parser.add_argument('-s', '--script_path', type=str, default='/data/jake/genefusion/scripts/python/fileid2sample_type.py', help='Path to file id to sample type script')
parser.add_argument('-l', '--log_file', type=str, default='/data/jake/genefusion/logs/migrate_sample2.log', help='Path to log file')
args = parser.parse_args()

with open(args.log_file, 'a') as log:
    # write time
    log.write(f'Running migrate_sample2.py at {time.ctime()}\n')
    # write args
    log.write(f'Input dir: {args.input_dir}\n')
    log.write(f'Normal dir: {args.normal_dir}\n')
    log.write(f'Tumour dir: {args.tumour_dir}\n')
    log.write(f'Script path: {args.script_path}\n')
    log.write(f'Log file: {args.log_file}\n')

# only lookup a file ids sample type once
def lookup_sample_type(file_id,cache):
    """
    Lookup the sample type (tumour/normal) for a given file id.
    """
    if file_id in cache.keys():
        print('file id found in cache:', file_id)
        return cache[file_id], cache
    try:
        print('looking up sample type for file id:', file_id)
        sample = subprocess.check_output([args.script_path, '-f', file_id]).decode('utf-8').strip()
        cache[file_id] = sample
        return sample, cache
    except subprocess.CalledProcessError as e:
        print("Command failed with error:", e)
        sys.exit(1)

valid = ['tumour', 'normal']
cache = {}
for f in os.listdir(args.input_dir):
    fid = os.path.basename(f).split('.')[0]

    # call the script to get sample type
    sample, cache = lookup_sample_type(fid,cache)
    if sample not in valid:
        with open(args.log_file, 'a') as log:
            log.write(f'File ID: {fid} Sample type: {sample} not recognised\n')
        continue

    # copy
    if sample == 'tumour':
        shutil.copy(os.path.join(args.input_dir,f), os.path.join(args.tumour_dir, f))
    elif sample == 'normal':
        shutil.copy(os.path.join(args.input_dir,f), os.path.join(args.normal_dir, f))

