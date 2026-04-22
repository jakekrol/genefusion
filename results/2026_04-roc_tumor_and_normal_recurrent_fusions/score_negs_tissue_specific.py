#!/usr/bin/env python


# positives are recurrent fusions from pcawg where
# we only use the evidence from tissues reported as recurrently present in pcawg
# to avoid the weighted sum from other tissues dominating the score
# example:
# fusion: erg-tmprss2
# tissue: prostate
# evidence types: prostate tumor, prostate normal, and 1000G

# negatives are recurrent fusions from pcawg where
# we use the evidence from a tissue that is not reported as recurrently present in pcawg
# example:
# fusion: erg-tmprss2
# tissue: breast
# evidence types: breast tumor, breast normal, and 1000G

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
parser.add_argument("--outfile", default="score_negs_tissue_specific.tsv")
parser.add_argument("--logfile", default="score_negs_tissue_specific.log")
parser.add_argument("--tumor_colmap", default="./tumor_colmap.yaml")
parser.add_argument("--normal_colmap", default="./normal_colmap.yaml")
parser.add_argument("--w_normals", default="0.25,0.5,0.75", type=str, help="csv string of w_normal weights to try")
parser.add_argument('--bed', default='../2025_04-gene_bedfile_cln/grch37.genes.bed')
args = parser.parse_args()

### inputs
df_fusions = pd.read_csv(args.fusions, sep="\t")
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
df_fusions["tissue_w_data"] = df_fusions['tissues'].apply(tissue_w_data)
# filter to fusions with at least one tissue with data
df_fusions = df_fusions[df_fusions["tissue_w_data"].apply(len) > 0]
df_fusions = df_fusions.reset_index(drop=True)
total= df_fusions.shape[0]
all_pcawg_tissues = df_evidence.columns.tolist()[2:]
# remove 1000g
all_pcawg_tissues = [t for t in all_pcawg_tissues if "1000g" not in t.lower()]
all_pcawg_tissues = set([t.split("_")[1].lower() for t in all_pcawg_tissues])
# join evidence onto fusions
df_fusions = left_sort_fusion_set(df_fusions, df_bed)
df_fusions = df_fusions.merge(df_evidence, on=['gene_left','gene_right'], how="left")
with open(args.logfile, "w") as f:
	f.write(f"{df_fusions.shape[0]} out of {total} fusions have at least one tissue with data\n")

def colmap_subset(colmap_tumor, colmap_normal, tissues, include_1000g=False):
	'''
	goal: select relevant evidence types for a fusion
	colmap_tumor: tumor evidence type dictionary
	colmap_normal: normal evidence type dictionary
	tissues: list of tissues the fusion is recurrent in 
	include_1000g: whether to include 1000g evidence
	'''
	colmap_t = {}
	for t in tissues:
		t = t.lower()
		for k, v in colmap_tumor.items():
			# get tissue of evidence cateogory
			tissue_key = k.split('_')[0].lower()
			if t == tissue_key:
				colmap_t[k] = v
	colmap_n = {}
	for t in tissues:
		t = t.lower()
		for k, v in colmap_normal.items():
			# get tissue of evidence cateogory
			tissue_key = k.split('_')[0].lower()
			if t == tissue_key:
				colmap_n[k] = v
	# optionally add 1000g evidence types
	if include_1000g:
		for k, v in colmap_normal.items():
			if "1000g" in k.lower():
				colmap_n[k] = v
	return colmap_t, colmap_n

### scoring positives
# for each fusion select only the relevance evidence types from the column maps
for i, row in df_fusions.iterrows():
	# try various w_normal weights for the score
	for w_normal in w_normals:
		tissues = row["tissue_w_data"]
		# evidence types for a fusion with respect to its recurrent tissues (for positive scoring)
		colmaps = []
		tumor_colmap_in_tissue, normal_colmap_in_tissue = colmap_subset(
			tumor_colmap,
			normal_colmap,
			tissues,
			include_1000g=False
		)
		colmaps.append((tumor_colmap_in_tissue, normal_colmap_in_tissue, f"score_as_positive-wnormal_{w_normal}"))
		tumor_colmap_in_tissue_w_1000g, normal_colmap_in_tissue_w_1000g = colmap_subset(
			tumor_colmap,
			normal_colmap,
			tissues,
			include_1000g=True
		)
		colmaps.append((tumor_colmap_in_tissue_w_1000g, normal_colmap_in_tissue_w_1000g, f"score_as_positive-w_1000g-wnormal_{w_normal}"))
		off_tissues = all_pcawg_tissues.difference(set(tissues))
		for t in off_tissues:
			# evidence types for using fusions as negative with respect to a specific tissue
			tumor_colmap_i_neg, normal_colmap_i_neg = colmap_subset(tumor_colmap, normal_colmap, [t], include_1000g=False)
			colmaps.append((tumor_colmap_i_neg, normal_colmap_i_neg, f"score_as_negative-tissue_{t}-wnormal_{w_normal}"))
			tumor_colmap_i_neg_w_1000g, normal_colmap_i_neg_w_1000g = colmap_subset(tumor_colmap, normal_colmap, [t], include_1000g=True)
			colmaps.append((tumor_colmap_i_neg_w_1000g, normal_colmap_i_neg_w_1000g, f"score_as_negative-tissue_{t}-w_1000g-wnormal_{w_normal}"))
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
	
     
     
     
