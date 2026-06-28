#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import pandas as pd
from polymerization.datasets import *
from polymerization.stix2fusion import *
from polymerization.io import *
import subprocess

parser = argparse.ArgumentParser(description='')
parser.add_argument('--score', default='score.in_tissue.tsv')
parser.add_argument('--outfile', default='hist_score.png')
parser.add_argument('--bed', default='../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed')
parser.add_argument("--density_script", default="density.py")
args = parser.parse_args()

df_bed = read_bed(args.bed, gene_col_idx=3)
df_score = pd.read_csv(args.score, sep='\t')
df_tumor_fusion = get_pcawg_recurrent_tumor_fusions()
df_tumor_fusion = left_sort_fusion_set(df_tumor_fusion, df_bed=df_bed)
df_score['gene_left'] = df_score['fusion'].apply(lambda x: x.split('--')[0])
df_score['gene_right'] = df_score['fusion'].apply(lambda x: x.split('--')[1])
df_merge = pd.merge(df_score, df_tumor_fusion, how='left', on=['gene_left', 'gene_right'])
df_merge_tumor = df_merge[df_merge['label']==1]
df_merge_tumor['reported_previously'] = df_merge_tumor['reported_previously'].apply(lambda x: 'novel' if x == 'False' else x)
ordering = ['chimerkb', 'chimerpub', 'chimerseq', 'novel']
files=[]
for cat in ordering:
	df_merge_tumor_cat = df_merge_tumor[df_merge_tumor['reported_previously']==cat]
	df_merge_tumor_cat['score'].to_csv(f"score.in_tissue.tumor_{cat}.tsv", sep='\t', index=False, header=False)
	print(f"Category: {cat}")
	print(df_merge_tumor_cat['score'].describe())
	files.append(f"score.in_tissue.tumor_{cat}.tsv")
df_merge_normal = df_merge[df_merge['label']==0]
df_merge_normal['score'].to_csv(f"score.in_tissue.normal.tsv", sep='\t', index=False, header=False)
files.append(f"score.in_tissue.normal.tsv")

subprocess.run(
    [args.density_script,
     "--inputs", ",".join(files), "--output", "density.in_tissue.png", "--names",
     "Tumor ChimerKB,Tumor ChimerPub,Tumor ChimerSeq,PCAWG novel,Normal fusions",
     "--title", "Recurrent tumor and normal fusion scores", "--ylabel", "Score"]
)


