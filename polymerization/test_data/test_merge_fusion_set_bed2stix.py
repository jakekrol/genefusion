#!/usr/bin/env python3
import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *
import time


df_fusionset=read_fusion_set('normal_fusion.modif.tsv')
df_bed=read_bed('grch37.genes.bed.added', gene_col_idx=3)
df_fusionset_sort = left_sort_fusion_set(df_fusionset, df_bed)
df_merged = merge_fusion_set_with_bed(df_fusionset_sort, df_bed)
df_shard=read_stix_shardfile('shardfile.tsv')

outdir='stix_output'
t_0=time.time()
cpus=30
merge_fusion_set_bed2stix(
	df_merged,
	df_shard,
	outdir,
	max_workers=cpus,
)
t_1=time.time()
print(f"Total time: {t_1-t_0} seconds")