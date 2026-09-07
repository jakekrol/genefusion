#!/usr/bin/env python3

import argparse
import pandas
from polymerization.datasets import *

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default="fusion_eval.tsv")
args = parser.parse_args()

COL_EVALUATION_DATASET="eval_dataset"

def merge(x,y, col_evaluation_dataset=COL_EVALUATION_DATASET):
    m = pd.merge(x,y, on = ['gene_left', 'gene_right'],how='outer')
    return m

def main():
    dataset='babiceanu_recurrent_normal_tissue_specific_fusions'
    df_babiceanu_recurrent_normal_tissue_specific_fusions = \
        get_babiceanu_recurrent_normal_tissue_specific_fusions()
    df_babiceanu_recurrent_normal_tissue_specific_fusions[COL_EVALUATION_DATASET] = dataset
    new_cols = ['gene_left', 'gene_right']
    for col in df_babiceanu_recurrent_normal_tissue_specific_fusions.columns:
        if not (col in ['gene_left', 'gene_right']):
            new_cols.append(f"{col}_{dataset}")
    df_babiceanu_recurrent_normal_tissue_specific_fusions.columns = new_cols

    dataset='babiceanu_recurrent_normal_tissue_agnostic_fusions'
    df_babiceanu_recurrent_normal_tissue_agnostic_fusions = \
        get_babiceanu_recurrent_normal_tissue_agnostic_fusions()
    df_babiceanu_recurrent_normal_tissue_agnostic_fusions[COL_EVALUATION_DATASET] = dataset
    new_cols = ['gene_left', 'gene_right']
    for col in df_babiceanu_recurrent_normal_tissue_agnostic_fusions.columns:
        if not (col in ['gene_left', 'gene_right']):
            new_cols.append(f"{col}_{dataset}")
    df_babiceanu_recurrent_normal_tissue_agnostic_fusions.columns = new_cols

    dataset='pcawg_recurrent_tumor_fusions'
    df_pcawg_recurrent_tumor_fusions = get_pcawg_recurrent_tumor_fusions()
    df_pcawg_recurrent_tumor_fusions[COL_EVALUATION_DATASET] = dataset
    new_cols = ['gene_left', 'gene_right']
    for col in df_pcawg_recurrent_tumor_fusions.columns:
        if not (col in ['gene_left', 'gene_right']):
            new_cols.append(f"{col}_{dataset}")
    df_pcawg_recurrent_tumor_fusions.columns = new_cols

    dataset='pcawg_tumor_fusions'
    df_pcawg_tumor_fusions = get_pcawg_tumor_fusions()
    df_pcawg_tumor_fusions[COL_EVALUATION_DATASET] = dataset
    new_cols = ['gene_left', 'gene_right']
    for col in df_pcawg_tumor_fusions.columns:
        if not (col in ['gene_left', 'gene_right']):
            new_cols.append(f"{col}_{dataset}")
    df_pcawg_tumor_fusions.columns = new_cols

    dataset='cosmic'
    df_cosmic_tumor_fusions = get_cosmic_tumor_fusions()
    df_cosmic_tumor_fusions[COL_EVALUATION_DATASET] = dataset
    new_cols = ['gene_left', 'gene_right']
    for col in df_cosmic_tumor_fusions.columns:
        if not (col in ['gene_left', 'gene_right']):
            new_cols.append(f"{col}_{dataset}")
    df_cosmic_tumor_fusions.columns = new_cols

    dataset='thousg_rna_star_fusion'
    df_thousg_rna_star_fusion = get_thousg_rna_star_fusion_calls()
    # consolidate all samples per fusion into single string
    thousg_rna_data = []
    for fusion, group in df_thousg_rna_star_fusion.groupby(["gene_left", "gene_right"]):
        samples = list(set(group['sample']))
        samples = ','.join(samples)
        thousg_rna_data.append((fusion[0], fusion[1], samples))
    df_thousg_rna_star_fusion = pd.DataFrame(thousg_rna_data, columns = ['gene_left','gene_right', 'samples'])
    df_thousg_rna_star_fusion[COL_EVALUATION_DATASET] = dataset
    new_cols = ['gene_left', 'gene_right']
    for col in df_thousg_rna_star_fusion.columns:
        if not (col in ['gene_left', 'gene_right']):
            new_cols.append(f"{col}_{dataset}")
    df_thousg_rna_star_fusion.columns = new_cols

    df_merge = merge(
        df_babiceanu_recurrent_normal_tissue_specific_fusions,
        df_babiceanu_recurrent_normal_tissue_agnostic_fusions
    )
    df_merge = merge(df_merge, df_pcawg_recurrent_tumor_fusions)
    df_merge = merge(df_merge, df_pcawg_tumor_fusions)
    df_merge = merge(df_merge, df_cosmic_tumor_fusions)
    df_merge = merge(df_merge, df_thousg_rna_star_fusion)
    df_merge.reset_index(drop=True, inplace=True)
    # consolidate eval data set columns
    eval_dataset_values = []
    for i, row in df_merge.iterrows():
        datasets=[]
        for col in df_merge.columns:
            if col.startswith('eval_dataset'):
                if not pd.isna(row[col]):
                    datasets.append(row[col])
        eval_dataset_values.append(
            ','.join(datasets)
        )
    # drop old eval dataset columns
    drop_cols = []
    for col in df_merge.columns:
        if col.startswith("eval_dataset"):
            drop_cols.append(col)
    df_merge.drop(columns=drop_cols, inplace=True)
    df_merge[COL_EVALUATION_DATASET] = eval_dataset_values
        
            
    leading_columns = ['gene_left', 'gene_right', COL_EVALUATION_DATASET]
    for col in df_merge.columns:
        if "tissue" in col:
            leading_columns.append(col)
    columns = df_merge.columns.tolist()
    columns = leading_columns + [c for c in columns if c not in leading_columns]
    df_merge = df_merge[columns]
    # df_merge.fillna('.', inplace=True)
    df_merge.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
    