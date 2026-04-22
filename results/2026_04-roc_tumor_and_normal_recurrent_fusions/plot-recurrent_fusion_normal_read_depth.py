#!/usr/bin/env python

import argparse
import os,sys
import pandas as pd
import subprocess
import ast
import tempfile


parser = argparse.ArgumentParser()
parser.add_argument("--input", '-i', default='./score_negs_tissue_specific.tsv')
parser.add_argument("--output", '-o', default='recurrent_fusion_normal_read_depth.png')
parser.add_argument("--density_script", default="/data/jake/rl-tools/plot/density.py")
parser.add_argument("--title", default="Normal read depth for recurrent fusions in/off tissue")
args = parser.parse_args()

df = pd.read_csv(parser.parse_args().input, sep="\t")
# for each row, lookup the normal evidence for the recurrent tissue
depths_in=[]
depths_off=[]
read_cols_normal= [ c for c in df.columns if c.startswith("reads_") and c.endswith("_normal_dna") ]
for i,row in df.iterrows():
    gene_left = row["gene_left"]
    gene_right = row["gene_right"]
    tissue = row['tissue_w_data']
    tissue = ast.literal_eval(tissue) # convert string representation of list to actual list
    # a fusion can have multiple recurrent tissues
    # most are just one
    for t in tissue:
        normal_col_in = f"reads_{t}_normal_dna"
        try:
            depth = row[normal_col_in]
            depths_in.append(depth)
        except KeyError:
            print(f"Warning: no normal read depth column for tissue {t} in row {i}")
            continue
        normal_cols_off = set(read_cols_normal) - set([normal_col_in])
        for normal_col in normal_cols_off:
            try:
                depth_off = row[normal_col]
                depths_off.append(depth_off)
            except KeyError:
                print(f"Warning: no normal read depth column {normal_col} in row {i}")
                continue
# write depths line delimited to tempfile
with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
    for d in depths_in:
        tmp.write(f"{d}\n")
    tmp_in = tmp.name
    print("# tmp: ", tmp_in)

with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
    for d in depths_off:
        tmp.write(f"{d}\n")
    tmp_off = tmp.name
    print("# tmp: ", tmp_off)

# plot histogram of depths
cmd = f"{args.density_script} -i {tmp_in},{tmp_off} -o {args.output} --ylabel 'Normal read depth' --title '{args.title}' --xlabel ' ' --names 'in,off' --show_median"
subprocess.run(cmd, shell=True)