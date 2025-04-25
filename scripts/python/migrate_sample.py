#!/usr/bin/env python3

import os,sys
import argparse
import subprocess
import shutil

#! this approach is really inefficient because we look up the sample over 20k times!
# instead let's just loop over files and cache the sample type of each file id

# input: 1) a sample unswapped intersect file, 2) normal dir, 3) tumour dir, 4) path to file id to sample type script
# output: None 
# effect: copy file to normal/tumour directories

parser = argparse.ArgumentParser(description='Copy sample intersect file to normal/tumour directories')
parser.add_argument('-i', '--intersect_file', type=str, required=True)
parser.add_argument('-n', '--normal_dir', type=str, required=True)
parser.add_argument('-t', '--tumour_dir', type=str, required=True)
parser.add_argument('-s', '--script_path', type=str, default='/data/jake/genefusion/scripts/python/fileid2sample_type.py', help='Path to file id to sample type script')
args = parser.parse_args()

# get file id
fid = os.path.basename(args.intersect_file).split('.')[0]
print('File ID:', fid)

# call the script to get sample type
try:
    sample = subprocess.check_output([args.script_path, '-f', fid]).decode('utf-8').strip()
    print('Sample type:', sample)
except subprocess.CalledProcessError as e:
    print("Command failed with error:", e)
assert sample in ['tumour', 'normal'], f'Sample type {sample} not recognised'

# copy
if sample == 'tumour':
    shutil.copy(args.intersect_file, os.path.join(args.tumour_dir, os.path.basename(args.intersect_file)))
if sample == 'normal':
    shutil.copy(args.intersect_file, os.path.join(args.normal_dir, os.path.basename(args.intersect_file)))

