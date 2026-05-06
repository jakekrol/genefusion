#!/usr/bin/env python3
import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *
from polymerization.giggle2fusion import *
import time

df_fusionset=read_fusion_set('normal_fusion.modif.tsv')
df_bed=read_bed('grch37.genes.bed.added', gene_col_idx=3)
df_fusionset_sort = left_sort_fusion_set(df_fusionset, df_bed)
df_merged = merge_fusion_set_with_bed(df_fusionset_sort, df_bed)
df_shard=read_giggle_shardfile('shardfile.giggle.tsv')
r1 = df_merged.iloc[0]
l = r1['gene_left']
r = r1['gene_right']
chrom_left = r1['chromosome_left']
start_left = r1['start_left']
end_left = r1['end_left']
index = df_shard.loc[0, 'giggle_index']
outfile = os.path.join(os.getcwd(), 'giggle_output.tsv')
argstring = f"search -i {index} -r {chrom_left}:{start_left}-{end_left} -v"
print(f"Running giggle with arguments: {argstring}")
print(f"Output will be written to: {outfile}")
run_giggle(argstring, l, r, outfile)

