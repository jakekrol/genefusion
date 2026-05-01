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
parser.add_argument("--fusions", default="./recurrent_tumor_fusions.tsv", help="Fusion set file")
parser.add_argument("--evidence", default="./all_categories-fusion_table.tsv", help="Evidence")
parser.add_argument("--shardfile", default="./shardfile.tsv", help="Shardfile")
parser.add_argument("--outfile", default="score.tsv")
parser.add_argument("--logfile", default="score.log")
parser.add_argument("--tumor_colmap", default="./tumor_colmap.yaml")
parser.add_argument("--normal_colmap", default="./normal_colmap.yaml")
# parser.add_argument("--include_1000g", action='store_true')
parser.add_argument("--w_normal", default=0.5, type=float, help="Weight for normal evidence in score calculation")
args = parser.parse_args()

### inputs
df_fusions = pd.read_csv(args.fusions, sep="\t")
df_evidence = pd.read_csv(args.evidence, sep="\t")
df_shardfile = read_stix_shardfile(args.shardfile)
with open(args.tumor_colmap) as f:
	tumor_colmap = yaml.safe_load(f)
with open(args.normal_colmap) as f:
	normal_colmap = yaml.safe_load(f)

### setup
# get evidence categories (e.g., tissues) from shardfile
categories = set(df_shardfile["category"])
# find which fusions we have data for their recurrent tissues
def tissue_w_data(tissues):
	t_w_data = []
	for t in tissues.split(","):
		for c in categories:
			# string check
			if t in c:
				t_w_data.append(t)
				break
	return t_w_data
df_fusions["tissue_w_data"] = df_fusions['tissues'].apply(tissue_w_data)
# filter to fusions with at least one tissue with data
total= df_fusions.shape[0]
df_fusions = df_fusions[df_fusions["tissue_w_data"].apply(len) > 0]
with open(args.logfile, "w") as f:
	f.write(f"{df_fusions.shape[0]} out of {total} fusions have at least one tissue with data\n")

### scoring
# for each fusion select only the relevance evidence types from the column maps
for i, row in df_fusions.iterrows():
	tissues = row["tissue_w_data"]
	# get tissue-specific tumor colmap for each fusion

	### colmaps as positive
	tumor_colmap_i_pos = {}
	for t in tissues:
		t = t.lower()
		for k,v in tumor_colmap.items():
			# get tissue of evidence cateogory
			tissue_key = k.split('_')[0].lower()
			if t == tissue_key:
				tumor_colmap_i_pos[k] = v
	# get tissue-specific normal colmap for each fusion
	normal_colmap_i_pos = {}
	for t in tissues:
		t = t.lower()
		for k,v in normal_colmap.items():
			# get tissue of evidence cateogory
			tissue_key = k.split('_')[0].lower()
			if t == tissue_key:
				normal_colmap_i_pos[k] = v
	
	normal_colmap_i_pos_w_1000g= normal_colmap_i_pos.copy()
	# add 1000G evidence types to normal colmap
	for k,v in normal_colmap.items():
		if ("high_cov" in k.lower()) and ("1000g" in k.lower()):
			normal_colmap_i_pos_w_1000g[k] = v
	### colmaps as negative
	tumor_colmap_i_neg = {}
	for t in tissues:
		t = t.lower()
		for k,v in tumor_colmap.items():
			tissue_key = k.split('_')[0].lower()
			# for negative fusions, we want to exclude the evidence from the recurrent tissue
			if t != tissue_key:
				tumor_colmap_i_neg[k] = v
	# as negative
	normal_colmap_i_neg = {}
	for t in tissues:
		t = t.lower()
		for k,v in normal_colmap.items():
			tissue_key = k.split('_')[0].lower()
			# no 1000g here
			if "1000g" in tissue_key:
				next
			if t != tissue_key:
				normal_colmap_i_neg[k] = v
	# as negative with 1000g
	normal_colmap_i_neg_w_1000g = normal_colmap_i_neg.copy()
	# add 1000G evidence types to normal colmap
	for k,v in normal_colmap.items():
		if ("high_cov" in k.lower()) and ("1000g" in k.lower()):
			normal_colmap_i_neg_w_1000g[k] = v
	### scoring
	## as positive
	# positive no 1000g
	T, N = fusion_tbl2_score_input(
		df_evidence.iloc[[i]],
		tumor_colmap_i_pos,
		normal_colmap_i_pos
	)
	y = score_numba_batched(T,N, w_normal=args.w_normal)[0]
	df_fusions.at[i, "score_as_positive"] = y
	# positive w 1000g
	T, N = fusion_tbl2_score_input(
		df_evidence.iloc[[i]],
		tumor_colmap_i_pos,
		normal_colmap_i_pos_w_1000g
	)
	y= score_numba_batched(T,N, w_normal=args.w_normal)[0]
	df_fusions.at[i, "score_as_positive_w_1000g"] = y
	# negative no 1000g
	T, N = fusion_tbl2_score_input(
		df_evidence.iloc[[i]],
		tumor_colmap_i_neg,
		normal_colmap_i_neg
	)
	y = score_numba_batched(T,N, w_normal=args.w_normal)[0]
	df_fusions.at[i, "score_as_negative"] = y
	# negative w 1000g
	T, N = fusion_tbl2_score_input(
		df_evidence.iloc[[i]],
		tumor_colmap_i_neg,
		normal_colmap_i_neg_w_1000g
	)
	y = score_numba_batched(T,N, w_normal=args.w_normal)[0]
	df_fusions.at[i, "score_as_negative_w_1000g"] = y

df_fusions = df_fusions.sort_values("score_as_positive", ascending=False)
df_fusions.to_csv(args.outfile, sep="\t", index=False)
	
     
     
     
