#!/usr/bin/env python3

from polymerization.io import *
from polymerization.analysis import *
import argparse
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Quantify inter and intra chromosomal fusion supporting reads')
parser.add_argument('--bed', default='../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed', help='Path to gene bed file')
parser.add_argument('--g2f_outdir', default='../2026_06-g2f-all_gene_pairs/g2f_out', help='Path to g2f output directory')
parser.add_argument('--outdir', default='.', help='Path to output directory')
parser.add_argument('--seed', type=int, default=0, help='Random seed for reproducibility')
parser.add_argument('--num_genes', type=int, default=1000, help='Number of genes to sample for analysis')
args = parser.parse_args()

# df_bed = read_bed(args.bed, gene_col_idx=3)
# dirs = os.listdir(args.g2f_outdir)
# dir_1000g = [d for d in dirs if 'high_coverage_1000g_dna' in d][0]
# dirs_normal = [d for d in dirs if 'normal' in d]
# dirs_tumor = [d for d in dirs if 'tumor' in d]
# genes = set(df_bed['gene_name'].sample(n=args.num_genes, random_state=args.seed).tolist())

# counts = {'gene': [], 'category': [], 'intra': [], 'inter': []}
# for i,gene in enumerate(genes):
# 	print(f"# gene {i+1}/{len(genes)}: {gene}")
# 	### 1000g
# 	inter_1000g=0
# 	intra_1000g=0
# 	path_1000g = os.path.join(args.g2f_outdir, dir_1000g, f"{gene}.giggle.clean.swap.intersect.bed.gz")
# 	if os.path.exists(path_1000g):
# 		try: 
# 			_, df_1000g = read_g2f_intersect(path_1000g, bgzip=True)
# 			data = intersect2chromosome_count_1000g = intersect2chromosome_count(df_1000g)
# 			inter_1000g += data['inter']
# 			intra_1000g += data['intra']
# 		except pd.errors.EmptyDataError:
# 			pass
# 	### normal
# 	inter_normal=0
# 	intra_normal=0
# 	for dir_normal in dirs_normal:
# 		path_normal = os.path.join(args.g2f_outdir, dir_normal, f"{gene}.giggle.clean.swap.intersect.bed.gz")
# 		if os.path.exists(path_normal):
# 			try:
# 				_, df_intersect_normal = read_g2f_intersect(path_normal, bgzip=True)
# 				data = intersect2chromosome_count(df_intersect_normal)
# 				inter_normal += data['inter']	
# 				intra_normal += data['intra']
# 			except pd.errors.EmptyDataError:
# 				pass
# 	### tumor
# 	inter_tumor=0
# 	intra_tumor=0
# 	for dir_tumor in dirs_tumor:
# 		path_tumor = os.path.join(args.g2f_outdir, dir_tumor, f"{gene}.giggle.clean.swap.intersect.bed.gz")
# 		if os.path.exists(path_tumor):
# 			try:
# 				_, df_intersect_tumor = read_g2f_intersect(path_tumor, bgzip=True)
# 				data = intersect2chromosome_count(df_intersect_tumor)
# 				inter_tumor += data['inter']	
# 				intra_tumor += data['intra']
# 			except pd.errors.EmptyDataError:
# 				pass
# 	counts['gene'].append(gene)
# 	counts['category'].append('1000g')
# 	counts['intra'].append(intra_1000g)
# 	counts['inter'].append(inter_1000g)
# 	counts['gene'].append(gene)
# 	counts['category'].append('normal')
# 	counts['intra'].append(intra_normal)
# 	counts['inter'].append(inter_normal)
# 	counts['gene'].append(gene)
# 	counts['category'].append('tumor')
# 	counts['intra'].append(intra_tumor)
# 	counts['inter'].append(inter_tumor)

# df_counts = pd.DataFrame(counts)
df_counts = pd.read_csv(os.path.join(args.outdir, 'chromosome_counts.tsv'), sep='\t')
df_counts['log_ratio_inter_div_intra'] = np.log10((df_counts['inter'] + 1) / (df_counts['intra'] + 1))
df_counts.to_csv(os.path.join(args.outdir, 'chromosome_counts.tsv'), index=False, sep='\t')

xmin=-5
xmax=2
fig, ax = plt.subplots(1,3, figsize=(10,5))
for i, category in enumerate(['1000g', 'normal', 'tumor']):
	subset = df_counts[df_counts['category'] == category]
	ax[i].hist(subset['log_ratio_inter_div_intra'], bins=20, color='black')
	ax[i].set_title(category)
	ax[i].set_xlabel('log10((inter + 1) / (intra + 1))', fontsize = 8)
	ax[i].set_ylabel('Num. genes')
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, 'log_ratio_histogram.category.png'))

# do a single histogram where we count gene-wise, not stratified by category
df_counts_grouped = df_counts.groupby('gene').agg({'inter': 'sum', 'intra': 'sum'}).reset_index()
plt.figure(figsize=(5,5))
df_counts_grouped['log_ratio_inter_div_intra'] = np.log10((df_counts_grouped['inter'] + 1) / (df_counts_grouped['intra'] + 1))
plt.hist(df_counts_grouped['log_ratio_inter_div_intra'], bins=20, color='black')
plt.title('All categories combined')
plt.xlabel('log10((inter + 1) / (intra + 1))', fontsize = 8)
plt.ylabel('Num. genes')
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, 'log_ratio_histogram.png'))
