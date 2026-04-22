#!/usr/bin/env python

import argparse
import subprocess
import pandas as pd
import tempfile
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", '-i', default='roc_data-neg_tissue_agnostic/roc_input-neg1-wnormal_0.0.tsv')
parser.add_argument("--output", '-o', default='set_sizes-neg_tissue_agnostic.bar.png')
parser.add_argument("--bar_script", default="/data/jake/rl-tools/plot/bars.py")
parser.add_argument("--title", default="Set sizes for tissue-agnostic experiment")
args= parser.parse_args()

df = pd.read_csv(args.input, sep="\t")
df['label'] = df['label'].replace({
	0: "Negative",
	1: "Positive"
})
with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
	for label, group in df.groupby("label"):
		size = group.shape[0]
		tmp.write(f"{label}\t{size}\n")
	tmp_path = tmp.name
	print("# tmp: ", tmp_path)

cmd = f"cat {tmp_path} | " \
    f"{args.bar_script} -o {args.output} --title '{args.title}'"
subprocess.run(cmd, shell=True)

os.remove(tmp_path)