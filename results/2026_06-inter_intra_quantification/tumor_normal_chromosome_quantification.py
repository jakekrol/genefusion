#!/usr/bin/env python3
from polymerization.io import *
from polymerization.analysis import *
from polymerization.datasets import *
import argparse
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description='Quantify inter and intra chromosomal fusion supporting reads')
parser.add_argument('--g2f_outdir', default='../2026_06-g2f-all_gene_pairs/g2f_out', help='Path to g2f output directory')
parser.add_argument('--outdir', default='.', help='Path to output directory')
args = parser.parse_args()


# dirs = os.listdir(args.g2f_outdir)
# dir_1000g = [d for d in dirs if 'high_coverage_1000g_dna' in d][0]
# dirs_normal = [d for d in dirs if 'normal' in d]
# dirs_tumor = [d for d in dirs if 'tumor' in d]

# # get tumor and normal gene sets
# df_normal = get_recurrent_normal_tissue_specific_fusions()
# df_tumor = get_pcawg_recurrent_tumor_fusions()
# genes_tumor = set(df_tumor['gene_x']).union(set(df_tumor['gene_y']))
# genes_normal = set(df_normal['gene_left']).union(set(df_normal['gene_right']))
# overlap = genes_tumor.intersection(genes_normal)
# genes_tumor = genes_tumor - overlap
# genes_normal = genes_normal - overlap


# ### normal fusion genes
# counts = {'gene': [], 'category': [], 'intra': [], 'inter': [], 'fusion_set': []}
# for i,gene in enumerate(genes_normal):
# 	print(f"# gene {i+1}/{len(genes_normal)}: {gene}")
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
# 	counts['fusion_set'].append('normal')

# 	counts['gene'].append(gene)
# 	counts['category'].append('normal')
# 	counts['intra'].append(intra_normal)
# 	counts['inter'].append(inter_normal)
# 	counts['fusion_set'].append('normal')

# 	counts['gene'].append(gene)
# 	counts['category'].append('tumor')
# 	counts['intra'].append(intra_tumor)
# 	counts['inter'].append(inter_tumor)
# 	counts['fusion_set'].append('normal')

# ### tumor fusion genes
# for i,gene in enumerate(genes_tumor):
# 	print(f"# gene {i+1}/{len(genes_tumor)}: {gene}")
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
# 	counts['fusion_set'].append('tumor')

# 	counts['gene'].append(gene)
# 	counts['category'].append('normal')
# 	counts['intra'].append(intra_normal)
# 	counts['inter'].append(inter_normal)
# 	counts['fusion_set'].append('tumor')

# 	counts['gene'].append(gene)
# 	counts['category'].append('tumor')
# 	counts['intra'].append(intra_tumor)
# 	counts['inter'].append(inter_tumor)
# 	counts['fusion_set'].append('tumor')


# df_counts = pd.DataFrame(counts)
# df_counts['log_ratio_inter_div_intra'] = np.log10((df_counts['inter'] + 1) / (df_counts['intra'] + 1))
# df_counts.to_csv(os.path.join(args.outdir, 'chromosome_counts.tumor_normal.tsv'), index=False, sep='\t')
df_counts = pd.read_csv(os.path.join(args.outdir, 'chromosome_counts.tumor_normal.tsv'), sep='\t')
df_counts_normal = df_counts[df_counts['fusion_set'] == 'normal']
df_counts_tumor = df_counts[df_counts['fusion_set'] == 'tumor']

xmin=-5
xmax=3
fig, ax = plt.subplots(2,3, figsize=(10,5))
# first row is normal fusion set, second row is tumor fusion set
for i, category in enumerate(['1000g', 'normal', 'tumor']):
	subset = df_counts_normal[df_counts_normal['category'] == category]
	ax[0,i].hist(subset['log_ratio_inter_div_intra'], bins=20, color='black')
	ax[0,i].set_title(f'{category}, normal fusion genes')
	ax[0,i].set_xlabel('log10((inter + 1) / (intra + 1))', fontsize = 8)
	ax[0,i].set_ylabel('Num. genes')
	ax[0,i].set_xlim(xmin, xmax)
# second row is tumor fusion set
for i, category in enumerate(['1000g', 'normal', 'tumor']):
	subset = df_counts_tumor[df_counts_tumor['category'] == category]
	ax[1,i].hist(subset['log_ratio_inter_div_intra'], bins=20, color='black')
	ax[1,i].set_title(f'{category}, tumor fusion genes')
	ax[1,i].set_xlabel('log10((inter + 1) / (intra + 1))', fontsize = 8)
	ax[1,i].set_ylabel('Num. genes')
	ax[1,i].set_xlim(xmin, xmax)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, 'log_ratio_histogram.category.tumor_normal.png'))

df_counts_grouped_normal = df_counts_normal.groupby('gene').agg({'inter': 'sum', 'intra': 'sum'}).reset_index()
plt.figure(figsize=(5,5))
df_counts_grouped_normal['log_ratio_inter_div_intra'] = np.log10((df_counts_grouped_normal['inter'] + 1) / (df_counts_grouped_normal['intra'] + 1))
plt.hist(df_counts_grouped_normal['log_ratio_inter_div_intra'], bins=20, color='black')
plt.title('Normal fusion genes')
plt.xlabel('log10((inter + 1) / (intra + 1))', fontsize = 8)
plt.ylabel('Num. genes')
plt.xlim(xmin, xmax)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, 'log_ratio_histogram.normal.png'))

df_counts_grouped_tumor = df_counts_tumor.groupby('gene').agg({'inter': 'sum', 'intra': 'sum'}).reset_index()
plt.figure(figsize=(5,5))
df_counts_grouped_tumor['log_ratio_inter_div_intra'] = np.log10((df_counts_grouped_tumor['inter'] + 1) / (df_counts_grouped_tumor['intra'] + 1))
plt.hist(df_counts_grouped_tumor['log_ratio_inter_div_intra'], bins=20, color='black')
plt.title('Tumor fusion genes')
plt.xlabel('log10((inter + 1) / (intra + 1))', fontsize = 8)
plt.ylabel('Num. genes')
plt.xlim(xmin, xmax)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, 'log_ratio_histogram.tumor.png'))