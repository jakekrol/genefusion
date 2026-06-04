#!/usr/bin/env python3

import argparse
import pandas as pd
from polymerization.io import *
from polymerization.stix2fusion import *

parser = argparse.ArgumentParser(description='Map normal fusions to tissues')
parser.add_argument('--input', default='../2026_06-normal_fusions_tissue_specific_set/recurrent_normal_fusion_set.filtered.tsv')	
parser.add_argument('--shardfile', default='./shardfile.tsv')
parser.add_argument('--bed', default='./grch37.genes.added.bed')
parser.add_argument('--output', default='normal_fusion_tissue_map.tsv')
args = parser.parse_args()

def sub_tissue_name(tissue):
	sub = {
     'lymph_node': 'lymph',
     'gall bladder': 'gallbladder',
     'bone_marrow': 'bone'

	}
	if tissue in sub:
		return sub[tissue]
	else:
		return tissue	

    

df_fusion = pd.read_csv(args.input, sep='\t')
df_bed = read_bed(args.bed, gene_col_idx=3)
df_shard = read_giggle_shardfile(args.shardfile)
categories = set(df_shard["category"])
data = []
for key, group in df_fusion.groupby(['up_gene', 'dw_gene']):
	gene_left = key[0]
	gene_right = key[1]
	tissues = list(group['sample_source'].unique())
	tissues.sort()
	data.append((gene_left, gene_right, ",".join(tissues)))

# make output table
df_output = pd.DataFrame(data, columns=['gene_x', 'gene_y', 'tissues'])
df_output_sort = left_sort_fusion_set(df_output, df_bed)
df_output_sort['tissues'] = df_output_sort['tissues'].str.lower()
df_output_sort.to_csv(args.output, sep='\t', index=False)
