#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import os
import pandas as pd
import sys
import yaml

parser = argparse.ArgumentParser()
# parser.add_argument("--indirs", help="csv string of directories with burden files")
parser.add_argument("--burden_file_suffix", default=".burden.txt")
parser.add_argument("--output", default="burden.tsv")
args = parser.parse_args()

# test it
def burdenfile2data(path_burden_file):
	burden=0
	gene = None
	with open(path_burden_file, "r") as f:
		for line in f:
			if line.startswith("#"):
				gene=line.split("=")[1].strip()
			else:
				burden = line.strip()
	return gene, burden
     
    

def main():
	print(f"# burden_file_suffix={args.burden_file_suffix}")
	print(f"# output={args.output}")
	# most datasets live here
	dirs = os.listdir("../2026_06-g2f-all_gene_pairs/g2f_out")
	dirs = [os.path.join("../2026_06-g2f-all_gene_pairs/g2f_out", d) for d in dirs]
	# 1kg rna lives here
	dirs.append("../2026_07-1kg-rna/g2f_out/mage_short_read_1000g_rna")
	gene2burden = defaultdict(dict)
	ndirs = len(dirs)
	for i,directory in enumerate(dirs):
		print (f"# scanning dir {i+1}/{ndirs}")
		files_burden = glob.glob(
			os.path.join(directory, f"*{args.burden_file_suffix}")
		)
		for f in files_burden:
			gene, burden = burdenfile2data(f)
			if gene:
				dataset = os.path.basename(os.path.dirname(f))
				gene2burden[gene][dataset] = burden

	columns = ["gene"] + [os.path.basename(d) for d in dirs]
	idx_col = {}
	for i, col in enumerate(columns):
		idx_col[col] = i
	data = []
	ngenes = len(gene2burden.keys())
	# make gene by dataset burden table
	for i,(gene, burden_dict) in enumerate(gene2burden.items()):
		if i % 100 == 0:
			print(f"# gene {i+1}/{ngenes}")
		row = [0] * (len(dirs)+1)
		row[0] = gene
		for col, burden in burden_dict.items():
			j = idx_col[col]
			row[j] = burden
		data.append(tuple(row))
	df = pd.DataFrame(data, columns = columns)
	df.to_csv(args.output,sep="\t",index=False)


      
				

       

if __name__ == "__main__":
    main()