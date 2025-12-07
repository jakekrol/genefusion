#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import shutil

# do it after excord bc files on aws
parser = argparse.ArgumentParser(description='Map file names to file IDs for excord files.')
parser.add_argument('-f', '--file_lookup', type=str, default='/data/jake/genefusion/results/2025_05-pcawg_filename2id/filename2fileid.tsv', help='File with mapping of filenames to file IDs')
parser.add_argument('-i', '--input', type=str, required=True, help='File with line delimited list of files to rename')
parser.add_argument('--suffix', type=str, default='.excord.bed.gz', help='Suffix to append to file IDs for output files')
parser.add_argument('--outdir', type=str, required=True, help='Output directory')
args = parser.parse_args()


df = pd.read_csv(args.file_lookup, sep='\t')
df['File_Name'] = df['File_Name'].apply(lambda x: x.split('.')[1])  # remove leading PCAWG_ and trailing part


with open(args.input, 'r') as f:
    files = [line.strip() for line in f if line.strip()]
for f in files:
    f_base = os.path.basename(f)
    f_key = f_base.split('.')[1]  # remove leading PCAWG_ and trailing part
    mask = df['File_Name'] == f_key
    if mask.any():
        fid = df[mask]['File_ID'].values[0]
        print(f'{f}\t{fid + args.suffix}')
        shutil.copy(f, os.path.join(args.outdir, fid + args.suffix))
    else:
        raise ValueError(f"File {f} not found in lookup table.")
