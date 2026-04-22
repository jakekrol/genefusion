#!/usr/bin/env python


import argparse
import yaml
from polymerization.io import *
from polymerization.score import *
from polymerization.stix2fusion import *
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--fusions", default="../2026_04-s2f_pcawg_recurrent/recurrent_tumor_fusions.tsv", help="Fusion set file")
parser.add_argument("--evidence", default="../2026_04-s2f_pcawg_recurrent/all_categories-fusion_table.tsv", help="Evidence")
parser.add_argument("--shardfile", default="../2026_04-s2f_pcawg_recurrent/shardfile.tsv", help="Shardfile")
parser.add_argument("--outfile", default="score-neg_tissue_agnostic.tsv")
parser.add_argument("--logfile", default="score-neg_tissue_agnostic.log")
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
df_fusions = pd.read_csv(args.fusions, sep="\t")
df_fusions.rename(columns={"tissues": "recurrent_tissues"}, inplace=True)
df_evidence = pd.read_csv(args.evidence, sep="\t")
df_shardfile = read_stix_shardfile(args.shardfile)
df_bed = read_bed(args.bed, gene_col_idx=3)
with open(args.tumor_colmap) as f:
	tumor_colmap = yaml.safe_load(f)
with open(args.normal_colmap) as f:
	normal_colmap = yaml.safe_load(f)
w_normals = [float(w) for w in args.w_normals.split(",")]

### setup
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
df_fusions["recurrent_tissues_w_data"] = df_fusions['recurrent_tissues'].apply(tissue_w_data)
# filter to fusions with at least one tissue with data
total= df_fusions.shape[0]
df_fusions = df_fusions[df_fusions["recurrent_tissues_w_data"].apply(len) > 0]
df_fusions = df_fusions.reset_index(drop=True)
# join evidence onto fusions
df_fusions = left_sort_fusion_set(df_fusions, df_bed)
df_fusions = df_fusions.merge(df_evidence, on=['gene_left','gene_right'], how="left")
with open(args.logfile, "w") as f:
	f.write(f"{df_fusions.shape[0]} out of {total} fusions have at least one tissue with data\n")

def colmap_subset(colmap_tumor, colmap_normal, tissues, include_1000g=False, match_tissue=True):
	'''
	goal: select relevant evidence types for a fusion
	colmap_tumor: tumor evidence type dictionary
	colmap_normal: normal evidence type dictionary
	tissues: list of tissues the fusion is recurrent in 
	include_1000g: whether to include 1000g evidence
	match_tissue: if true, then include tissues where the fusion is recurrent (positive scoring); if false, then include tissues where the fusion is not recurrent (for negative scoring)
	'''
	colmap_t = {}
	for t in tissues:
		t = t.lower()
		for k, v in colmap_tumor.items():
			# get tissue of evidence cateogory
			tissue_key = k.split('_')[0].lower()
			if match_tissue:
				if t == tissue_key:
					colmap_t[k] = v
			elif not match_tissue:
				if t != tissue_key:
					colmap_t[k] = v
	colmap_n = {}
	for t in tissues:
		t = t.lower()
		for k, v in colmap_normal.items():
			# get tissue of evidence cateogory
			tissue_key = k.split('_')[0].lower()
			if match_tissue:
				if t == tissue_key:
					colmap_n[k] = v
			elif not match_tissue:
				if t != tissue_key:
					colmap_n[k] = v
	# optionally add 1000g evidence types
	if include_1000g:
		for k, v in colmap_normal.items():
			if "1000g" in k.lower():
				colmap_n[k] = v
	return colmap_t, colmap_n

### scoring
# for each fusion select only the relevance evidence types from the column maps
for i, row in df_fusions.iterrows():
	tissues = row["recurrent_tissues_w_data"]
	# evidence types for using fusions as either positive or negative 
	tumor_colmap_i_pos, normal_colmap_i_pos = colmap_subset(tumor_colmap, normal_colmap, tissues, include_1000g=False)
	tumor_colmap_i_pos_w_1000g, normal_colmap_i_pos_w_1000g = colmap_subset(tumor_colmap, normal_colmap, tissues, include_1000g=True)
	tumor_colmap_i_neg, normal_colmap_i_neg = colmap_subset(tumor_colmap, normal_colmap, tissues, include_1000g=False, match_tissue=False)
	tumor_colmap_i_neg_w_1000g, normal_colmap_i_neg_w_1000g = colmap_subset(tumor_colmap, normal_colmap, tissues, include_1000g=True, match_tissue=False)
	# try various w_normal weights for the score
	for w_normal in w_normals:
		colmaps = [
			(tumor_colmap_i_pos, normal_colmap_i_pos, f"score_as_positive-wnormal_{w_normal}"),
			(tumor_colmap_i_pos_w_1000g, normal_colmap_i_pos_w_1000g, f"score_as_positive-w_1000g-wnormal_{w_normal}"),
			(tumor_colmap_i_neg, normal_colmap_i_neg, f"score_as_negative-wnormal_{w_normal}"),
			(tumor_colmap_i_neg_w_1000g, normal_colmap_i_neg_w_1000g, f"score_as_negative-w_1000g-wnormal_{w_normal}")
		]
		# make input arrays for scorying
		for tumor_colmap_i, normal_colmap_i, score_col in colmaps:
			T, N = fusion_tbl2_score_input(
				# row of evidence for the fusion
				df_fusions.iloc[[i]],
				tumor_colmap_i,
				normal_colmap_i
			)
			y = score_numba_batched(T,N, w_normal=w_normal)[0]
			df_fusions.at[i, score_col] = y
df_fusions.to_csv(args.outfile, sep="\t", index=False)
	
     
     
     

