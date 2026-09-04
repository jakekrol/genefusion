#!/usr/bin/env python3

import argparse
import pandas
from polymerization.datasets import *

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default="fusion_eval.tsv")
args = parser.parse_args()

COL_EVALUATION_DATASET="eval_dataset"

def merge(x,y, suffixes):
    return pd.merge(x,y, on = ['gene_left', 'gene_right'],how='outer', suffixes=suffixes)

def main():
	df_babiceanu_recurrent_normal_tissue_specific_fusions = \
		get_babiceanu_recurrent_normal_tissue_specific_fusions()
	df_babiceanu_recurrent_normal_tissue_specific_fusions[COL_EVALUATION_DATASET] = \
		'babiceanu_recurrent_normal_tissue_specific_fusions'

	df_babiceanu_recurrent_normal_tissue_agnostic_fusions = \
		get_babiceanu_recurrent_normal_tissue_agnostic_fusions()
	df_babiceanu_recurrent_normal_tissue_agnostic_fusions[COL_EVALUATION_DATASET] = \
		'babiceanu_recurrent_normal_tissue_agnostic_fusions'

	df_pcawg_recurrent_tumor_fusions = get_pcawg_recurrent_tumor_fusions()
	df_pcawg_recurrent_tumor_fusions[COL_EVALUATION_DATASET] = 'pcawg_recurrent_tumor_fusions'

	df_pcawg_tumor_fusions = get_pcawg_tumor_fusions()
	df_pcawg_tumor_fusions[COL_EVALUATION_DATASET] = 'pcawg_tumor_fusions'

	df_cosmic_tumor_fusions = get_cosmic_tumor_fusions()
	df_cosmic_tumor_fusions[COL_EVALUATION_DATASET] = 'cosmic'

	df_merge = merge(
     df_babiceanu_recurrent_normal_tissue_specific_fusions,
     df_babiceanu_recurrent_normal_tissue_agnostic_fusions,
     suffixes=('_babiceanu_recurrent_normal_tissue_specific_fusions', '_babiceanu_recurrent_normal_tissue_agnostic_fusions'))
	df_merge = merge(df_merge, df_pcawg_recurrent_tumor_fusions, suffixes=('', '_pcawg_recurrent_tumor_fusions'))
	df_merge = merge(df_merge, df_pcawg_tumor_fusions, suffixes=('', '_pcawg_tumor_fusions'))
	df_merge = merge(df_merge, df_cosmic_tumor_fusions, suffixes=('', '_cosmic'))
	leading_columns = ['gene_left', 'gene_right', COL_EVALUATION_DATASET]
	columns = df_merge.columns.tolist()
	columns = leading_columns + [c for c in columns if c not in leading_columns]
	df_merge = df_merge[columns]
	# consolidate the COL_EVALUATION_DATASET column to a single value per row
	df_merge.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
    