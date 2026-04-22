#!/usr/bin/env python
import argparse
import subprocess
import pandas as pd
import os,sys
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--roc-script", default="/data/jake/rl-tools/plot/roc.py",
                    help="Path to roc plotting script")
parser.add_argument("--score-table", default="score_negs_tissue_specific.tsv", help="Path to score table tsv file")
parser.add_argument("--outdir", default="roc_negs_tissue_specific")
parser.add_argument("--w_normals", default="0.25,0.5,0.75", type=str, help="csv string of w_normal weights to try")
args = parser.parse_args()
w_normals = [float(w) for w in args.w_normals.split(",")]
os.makedirs(args.outdir, exist_ok=True)

df= pd.read_csv(args.score_table, sep="\t")
for w_normal in w_normals:
	for i in [0,1]:
		outfile = os.path.join(args.outdir, f"roc_input-neg2-w_normal_{w_normal}")
		if i == 0:
			thousg=True
			outfile += "-w_1000g.tsv"
		else:
			thousg=False
			outfile += ".tsv"
		with open(outfile, "w") as f:
			f.write("gene_left\tgene_right\ttissue\tscore\tlabel\tscore_column_name\n")
			for i,row in df.iterrows():
				# drop nan values
				# these are the negative score columns for the tissue that the fusion is recurrent in.
				mask = row.isna()
				row = row[~mask]
				gene_left=row["gene_left"]
				gene_right=row["gene_right"]
				tissue=row['tissue_w_data']
				p_idx = [idx for idx in row.index if 'score_as_positive' in idx]
				p_idx = [idx for idx in p_idx if str(w_normal) in idx]
				n_idx = [idx for idx in row.index if 'score_as_negative' in idx]
				n_idx = [idx for idx in n_idx if str(w_normal) in idx]
				if not thousg:
					n_idx = [idx for idx in n_idx if "1000g" not in idx]
					p_idx = [idx for idx in p_idx if "1000g" not in idx]
				for idx in p_idx:
					score= row[idx]
					label=1
					f.write(f"{gene_left}\t{gene_right}\t{tissue}\t{score}\t{label}\t{idx}\n")
				for idx in n_idx:
					score= row[idx]
					label=0
					f.write(f"{gene_left}\t{gene_right}\t{tissue}\t{score}\t{label}\t{idx}\n")

# plot roc curve
rocfiles = glob.glob(os.path.join(args.outdir, "roc_*.tsv"))
print(rocfiles)
rocfiles = [os.path.basename(f) for f in rocfiles]
rocfiles = ",".join(rocfiles)
cmd = f"cd {args.outdir} && python {args.roc_script} --scores {rocfiles} " \
    f"--output roc_curves.png --score_col 3 --label_col 4 --header"
subprocess.run(cmd, shell=True)
