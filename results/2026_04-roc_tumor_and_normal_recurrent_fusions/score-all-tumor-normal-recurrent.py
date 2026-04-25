#!/usr/bin/env python


import argparse
import yaml
from polymerization.io import *
from polymerization.score import *
from polymerization.stix2fusion import *
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--tumor-fusions", default="../2026_04-s2f_pcawg_recurrent/recurrent_tumor_fusions.tsv", help="Fusion set file")
parser.add_argument("--normal-fusions", default="../2026_04-s2f_babiceanu_recurrent_normal_tissue_agnostic/recurrent_normal_fusions.tsv", help="Fusion set file")
parser.add_argument("--tumor-evidence", default="../2026_04-s2f_pcawg_recurrent/all_categories-fusion_table.tsv", help="Evidence")
parser.add_argument("--normal-evidence", default="../2026_04-s2f_babiceanu_recurrent_normal_tissue_agnostic/all_categories-fusion_table.tsv", help="Evidence")
parser.add_argument("--shardfile", default="../2026_04-s2f_pcawg_recurrent/shardfile.tsv", help="Shardfile")
parser.add_argument("--outfile", default="score-all-tumor-normal-recurrent.tsv")
parser.add_argument("--logfile", default="score-all-tumor-normal-recurrent.log")
parser.add_argument("--tumor_colmap", default="./tumor_colmap.yaml")
parser.add_argument("--normal_colmap", default="./normal_colmap.yaml")
parser.add_argument("--w_normals", default="0,0.25,0.5,0.75,1", type=str, help="csv string of w_normal weights to try")
parser.add_argument('--bed', default='../2025_04-gene_bedfile_cln/grch37.genes.bed')
args = parser.parse_args()

# print all args and their values
print("# arguments:")
for arg in vars(args):
	print(f"# {arg}: {getattr(args, arg)}")

### inputs
df_shardfile = read_stix_shardfile(args.shardfile)
df_bed = read_bed(args.bed, gene_col_idx=3)
# tumor
df_fusions_tumor = pd.read_csv(args.tumor_fusions, sep="\t")
df_fusions_tumor.rename(columns={"tissues": "recurrent_tissues"}, inplace=True)
df_fusions_tumor['label'] = 1
df_evidence_tumor = pd.read_csv(args.tumor_evidence, sep="\t")
# normal
df_fusions_normal = read_fusion_set(args.normal_fusions, sep="\t")
df_fusions_normal.rename(columns={"tissues": "recurrent_tissues"}, inplace=True)
df_fusions_normal['label'] = 0
df_evidence_normal = pd.read_csv(args.normal_evidence, sep="\t")
# get evidence categories (e.g., tissues) from shardfile
categories = set(df_shardfile["category"])
# find which fusions we have data for their recurrent tissues
def tissue_w_data(tissues):
	t_w_data = []
	for t in tissues.split(","):
		for c in categories:
			category_tissue = c.split("_")[0].lower()
			# string check
			if t == category_tissue:
				t_w_data.append(t)
				break
	return t_w_data
df_fusions_tumor["recurrent_tissues_w_data"] = df_fusions_tumor['recurrent_tissues'].apply(tissue_w_data)
# filter to fusions with at least one tissue with data
total= df_fusions_tumor.shape[0]
df_fusions_tumor = df_fusions_tumor[df_fusions_tumor["recurrent_tissues_w_data"].apply(len) > 0]
# combine tumor and normal fusion sets and evidence
df_fusions = pd.concat([df_fusions_tumor, df_fusions_normal], ignore_index=True)
df_evidence = pd.concat([df_evidence_tumor, df_evidence_normal], ignore_index=True)
# merge evidence onto fusions
df_fusions = left_sort_fusion_set(df_fusions, df_bed)
df_fusions = pd.merge(df_fusions, df_evidence, on=['gene_left', 'gene_right'], how='left')
with open(args.tumor_colmap) as f:
	tumor_colmap = yaml.safe_load(f)
with open(args.normal_colmap) as f:
	normal_colmap = yaml.safe_load(f)
w_normals = [float(w) for w in args.w_normals.split(",")]

### scoring
for i, row in df_fusions.iterrows():
	# try various w_normal weights for the score
	for w_normal in w_normals:
		for j in [0,1]: # thousand genomes indicator
			if j == 0:
				# drop 1000G from the normal_colmap
				normal_colmap_j = {k:v for k,v in normal_colmap.items() if not ("1000g" in k)}
				score_col_name = f"score-wnormal_{w_normal}"
			else:
				normal_colmap_j = normal_colmap
				score_col_name = f"score-w_1000g-wnormal_{w_normal}"
			# create score column
			T,N = fusion_tbl2_score_input(
				df_fusions.iloc[[i]],
				tumor_colmap,
				normal_colmap_j
			)
			y = score_numba_batched(T,N, w_normal=w_normal)[0]
			df_fusions.at[i, score_col_name] = y
df_fusions.to_csv(args.outfile, sep="\t", index=False)
	
     
     
     

