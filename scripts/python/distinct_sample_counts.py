#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse
import time

# input: a giggle intersect unswapped file
# output: a 3 col tsv of gene_left, gene_right, and number of samples >= 1 PE evidence of fusion

parser = argparse.ArgumentParser(description='Count samples with gene fusions')
parser.add_argument('-i', '--input', type=str, required=True, help='Input file of gene-wise fusion count files')
parser.add_argument('-l', '--left', type=str,  help='Left gene name')
parser.add_argument('-s', '--sample_col_idx', type=int, default=15, help='Sample column index (1-indexed)')
parser.add_argument('-r', '--right_gene_col_idx', type=int, default=4, help='Right gene column index (1-indexed)')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for sample count data')
parser.add_argument('-z', '--hack', action='store_true', help='Hack to parse left gene name from the input file name')
args = parser.parse_args()

args.sample_col_idx -= 1
args.right_gene_col_idx -= 1

if args.hack:
    args.left = '.'.join(os.path.basename(args.input).split('.')[:-4]) # remove the last 4 elements (chrom, strand, start, end)
    assert len(args.left) > 0, f'Left gene name could not be parsed from the input file name: result={args.left}. Please provide a left gene name using -l option.'
    print(f'Parsed left gene name from input file name: {args.left}')

t = time.time()
df = pd.read_csv(args.input, sep='\t', header=None, usecols=[args.right_gene_col_idx, args.sample_col_idx])
df.columns = ['right', 'sample_count']
x=df.groupby('right')['sample_count'].nunique()
x = x.to_frame()
x.reset_index(inplace=True) # make right a column
x['left'] = args.left
x = x[['left', 'right', 'sample_count']]
x.to_csv(args.output, sep='\t', header=True,index=False)
print(f'Finished counting distinct samples for fusions with left={args.left} in {args.input} in {time.time()-t:.2f} seconds')
