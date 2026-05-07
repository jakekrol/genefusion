#!/usr/bin/env python3
import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *
from polymerization.giggle2fusion import *
import time

# OUTDIR='g2f_output'
# df_fusionset=read_fusion_set('normal_fusion.modif.tsv')
# df_bed=read_bed('grch37.genes.bed.added', gene_col_idx=3)
# df_fusionset_sort = left_sort_fusion_set(df_fusionset, df_bed)
# df_merged = merge_fusion_set_with_bed(df_fusionset_sort, df_bed)
# df_shard=read_giggle_shardfile('shardfile.giggle.tsv')
# # giggle queries
# df_giggle = merge_fusion_set_bed2giggle(
#     df_merged, df_shard, outdir=OUTDIR, max_workers=30, timeout=60*60*2, bgzip=True
# )
# print(df_giggle.head())

# # swap intervals in giggle output
# df_swap = giggle2swap(df_giggle, bgzip=True)
# print(df_swap.head())

# bedtools intersect
path_bed='grch37.genes.bed.added'
path_giggle_swap='g2f_output/low_coverage_1000g/TMEM120A.giggle.swap.gz'
outfile='TMEM120A.giggle.swap.intersect.bed.gz'
bedtools_bin='/data/jake/bedtools.static.binary'
bedtools_intersect(path_giggle_swap, path_bed, outfile, bgzip=True, bedtools_bin=bedtools_bin)

