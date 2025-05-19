#!/usr/bin/env python3

import pandas as pd
import os
import sys
import argparse
import time

# input: a dir of w/ specimen stratified intersect files
# output: a 3 col tsv of gene_left, gene_right, and number of samples with the fusion

SAMPLECOLIDX=14
RIGHTGENECOLIDX=3

parser = argparse.ArgumentParser(description='Count samples with gene fusions')
parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input directory of gene-wise fusion count files')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for gene fusion table')
parser.add_argument('-l', '--logfile', type=str, default='count_samples_w1.log', help='Log file for gene fusion table')
args = parser.parse_args()

if os.path.exists(args.output):
    print(f'Output file {args.output} already exists. Exiting.')
    sys.exit(1)
files = [ os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir)]
assert len(files) > 0, f'No files found in {args.input_dir}'

with open(args.logfile, 'w') as log:
    # print the time
    log.write(f'Log file created at {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
    log.write(f'Input directory: {args.input_dir}\n')
    log.write(f'Output file: {args.output}\n')
n = len(files)
i = 1
for f in files:
    if os.path.getsize(f) == 0:
        print(f'File {f} is empty. Skipping.')
        with open(args.logfile, 'a') as log:
            log.write(f'File {f} is empty. Skipping.\n')
        continue
    left = os.path.basename(f).split('.')[:-4] # remove the last 4 elements (chrom, strand, start, end)
    left = '.'.join(left)
    df_f = pd.read_csv(f, sep='\t', header=None, usecols=[RIGHTGENECOLIDX, SAMPLECOLIDX])
    df_r_sample_count = df_f.groupby(RIGHTGENECOLIDX)[SAMPLECOLIDX].apply(lambda x: len(set(x)))
    for right, sample_count in df_r_sample_count.items():
        with open(args.output, 'a') as out:
            out.write(f'{left}\t{right}\t{sample_count}\n')
    print(f'Processed {i}/{n} files')
    i += 1



