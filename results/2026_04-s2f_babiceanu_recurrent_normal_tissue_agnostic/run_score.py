#!/usr/bin/env python
from polymerization.score import *
from polymerization.stix2fusion import *
import pandas as pd
import time
import yaml
UPPER_BOUND = 100
df_fusion = pd.read_csv('all_categories-fusion_table.tsv', sep='\t')
with open("tumor_colmap.yaml", 'r') as f:
	tumor_colmap = yaml.safe_load(f)
with open("normal_colmap.yaml", 'r') as f:
	normal_colmap = yaml.safe_load(f)
T, N = fusion_tbl2_score_input(
	df_fusion,
	tumor_colmap,
	normal_colmap
)
fusions=df_fusion['gene_left'] + '--' + df_fusion['gene_right']
t_0=time.time()
y = score_numba_batched(T,N)
t_1 = time.time()
print(f"score_numba_batched took {t_1-t_0:.2f} seconds")
for i,y_i in enumerate(y):
	fusion = fusions[i]
	score = y_i
	print(f"{fusion}\t{score:.4f}")