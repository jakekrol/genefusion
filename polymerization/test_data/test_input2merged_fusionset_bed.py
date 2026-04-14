#!/usr/bin/env python

import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *


df_fusionset=read_fusion_set('normal_fusion.modif.tsv')
df_bed=read_bed('grch37.genes.bed.added', gene_col_idx=3)
df_fusionset_sort = left_sort_fusion_set(df_fusionset, df_bed)
df_merged = merge_fusion_set_with_bed(df_fusionset_sort, df_bed)
print(df_merged.head())
