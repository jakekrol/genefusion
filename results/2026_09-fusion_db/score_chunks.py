#!/usr/bin/env python3
import numpy as np
import pandas as pd
from polymerization.score import coverage_normalize_df_evidence, coverage_normalize_tumor_reads, coverage_normalize_normal_reads, normalize_samples, burden_normalize_df_evidence, burden_normalize_reads
from polymerization.datasets import get_pcawg_data_types
import sqlite3
import time
import yaml

DATABASE="fusion.db"
TABLE='fusion_evidence'
TABLE_BURDEN='burden'
COLUMN_MAP='score_column_map.yaml'
NORMAL_COLUMN_PATTERNS = ["normal", "1000g", "thousg", "mage"]
THOUSG_COLUMN_PATTERNS = ["1000g", "thousg", "mage"]
with open(COLUMN_MAP, 'r') as f:
    column_map = yaml.safe_load(f)


def column_map2subpopulation_weights(column_map: dict):
    total_samples_normal=0
    total_samples_tumor=0
    subpopulation_weights = dict()
    # first get total tumor/normal samples
    for column,dict_metadata in column_map.items():
        evidence_type = dict_metadata['evidence_type']
        if evidence_type == 'sample':
            specimen = dict_metadata['specimen']
            samples = dict_metadata['total_samples']
            if specimen == 'normal':
                total_samples_normal += samples
            if specimen == 'tumor':
                total_samples_tumor += samples
    # weigh each subpopulation by fractional sample contribution to tumor/normal
    for column,dict_metadata in column_map.items():
        specimen = dict_metadata['specimen']
        samples = dict_metadata['total_samples']
        if specimen == 'normal':
            subpopulation_weights[column] = samples / total_samples_normal
        if specimen == 'tumor':
            subpopulation_weights[column] = samples / total_samples_tumor
    return subpopulation_weights

def get_pcawg_tumor_tissues():
    data = get_pcawg_data_types()
    tissues = set()
    for tissue, _ in data.items():
        if 'tumor' in data[tissue]['modality']['dna']['specimen']:
            tissues.add(tissue)
        if 'tumor' in data[tissue]['modality']['rna']['specimen']:
            tissues.add(tissue)
    return tissues

def chunk2df(
    database: str,
    table: str,
    table_burden: str,
    subpopulation_weights: dict,
    column_map: dict,
    tumor_tissues: set,
    chunk_size: int = int(10**5)
):
    # burden
    conn = sqlite3.connect(database)
    df_burden = pd.read_sql_query(f"SELECT * FROM {table_burden}", conn)
    df_burden.set_index("gene", inplace=True)
    n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    k_chunks = n // chunk_size # quotient
    if n % chunk_size != 0:
        k_chunks += 1 # add 1 if there is remainder
    times = []
    t_0 = time.time()
    # to-do
    # 1. wrap for loop in "with conn:" block to execute single transaction
    # 2. write directly to sql, not text files on disk
    with conn:
        for i_chunk,df in enumerate(pd.read_sql_query(f"SELECT * FROM {table}", conn, chunksize=chunk_size)):
            print(f"# scoring chunk {i_chunk+1}/{k_chunks}")
            # normalize
            df_norm_coverage = coverage_normalize_df_evidence(df.copy(), column_map) # func operates in place, use a copy here, bc still need input for next function
            df_norm_burden = burden_normalize_df_evidence(df, df_burden, column_map)
            ### scale
            # each read/sample normalized column gets a 1/2 scaling. So, when read+sample is summed it's in [0,1]
            # when combining across subpopulations, a weighted average is used. the summed result is also in [0,1]
            # if you write out the equation, each column is just multiplied by f/2. where f in [0,1] is the subpopulation weight
            # final score is adding the normalized and scaled columns
            ## global score
            # just apply weights and scaling to all columns, then add
            df_norm_coverage_global = df_norm_coverage.copy() # copy here, bc we need another copy for tissue-wise scoring
            df_norm_burden_global = df_norm_burden.copy()
            for column,subpop_w in subpopulation_weights.items():
                specimen = column_map[column]['specimen']
                if specimen == 'normal':
                    sign = -1
                elif specimen == 'tumor':
                    sign = 1
                else:
                    raise ValueError(f"specimen={specimen} for column={column} must be 'tumor' or 'normal'")
                df_norm_coverage_global[column] = df_norm_coverage_global[column] * subpop_w * 0.5 * sign
                df_norm_burden_global[column] = df_norm_burden_global[column] * subpop_w * 0.5 * sign
            df_norm_coverage_global['score_global_norm_coverage'] = df_norm_coverage_global[subpopulation_weights.keys()].sum(axis=1)
            df_norm_burden_global['score_global_norm_burden'] = df_norm_burden_global[subpopulation_weights.keys()].sum(axis=1)
            # horizontally stack all score columns
            df_out = df_norm_coverage_global[['gene_left', 'gene_right']]
            for df in [df_norm_coverage_global, df_norm_burden_global]:
                # select score columns
                columns = []
                for column in df.columns:
                    if 'score' in column:
                        assert (df[column].min() >= -1) and (df[column].max() <= 1)
                        columns.append(column)
                df_out = pd.concat([df_out, df[columns]],axis=1)
            ## tissue-wise scoring
            # need to adjust subpopulation weights here because a subset of populations (tissue specific plus 1kg) are used
            for tissue in tumor_tissues:
                tumor_cols = []
                normal_cols = []
                for column in column_map.keys():
                    specimen = column_map[column]['specimen']
                    # tissue specific
                    if (tissue in column) and specimen == 'tumor':
                        tumor_cols.append(column)
                    if (tissue in column) and specimen == 'normal':
                        normal_cols.append(column)
                    # thousg (control)
                    for pattern in THOUSG_COLUMN_PATTERNS:
                        if pattern in column:
                            normal_cols.append(column)
                # get adjusted subpopulation weights
                column_map_subset = dict()
                for k,v in column_map.items():
                    if (k in tumor_cols) or (k in normal_cols):
                        column_map_subset[k] = v
                subpopulation_weights_tissue = column_map2subpopulation_weights(column_map_subset)
                # subset relevant evidence columns and fusion keys
                df_norm_coverage_tissue = df_norm_coverage[['gene_left', 'gene_right'] + list(subpopulation_weights_tissue.keys())]
                df_norm_burden_tissue = df_norm_burden[['gene_left', 'gene_right'] + list(subpopulation_weights_tissue.keys())]
                for column,subpop_w in subpopulation_weights_tissue.items():
                    specimen = column_map[column]['specimen']
                    if specimen == 'normal':
                        sign = -1
                    elif specimen == 'tumor':
                        sign = 1
                    else:
                        raise ValueError(f"specimen={specimen} for column={column} must be 'tumor' or 'normal'")
                    df_norm_coverage_tissue[column] = df_norm_coverage_tissue[column] * subpop_w * 0.5 * sign
                    df_norm_burden_tissue[column] = df_norm_burden_tissue[column] * subpop_w * 0.5 * sign
                df_norm_coverage_tissue[f'score_{tissue}_norm_coverage'] = df_norm_coverage_tissue[subpopulation_weights_tissue.keys()].sum(axis=1)
                df_norm_burden_tissue[f'score_{tissue}_norm_burden'] = df_norm_burden_tissue[subpopulation_weights_tissue.keys()].sum(axis=1)
                # horizontally stack all score columns
                for df in [df_norm_coverage_tissue, df_norm_burden_tissue]:
                    columns = []
                    for column in df.columns:
                        if 'score' in column:
                            assert (df[column].min() >= -1) and (df[column].max() <= 1)
                            columns.append(column)
                    df_out = pd.concat([df_out, df[columns]],axis=1)
            if_exists = "replace" if i_chunk == 0 else "append"
            df_out.to_sql(
                "score",
                conn,
                if_exists=if_exists,
                index=False,
            )
            t_elapsed = time.time() - t_0
            print(f"# time for chunk {i_chunk+1}: {t_elapsed}")
            times.append(t_elapsed)
            t_0 = time.time()
        conn.close()
    return np.mean(times), chunk_size, n



# def estimate_runtime(avg_t_per_chunk, chunk_size, n):
#     iterations = n / chunk_size
#     t_est = avg_t_per_chunk * iterations
#     return t_est

def main():
    subpopulation_weights = column_map2subpopulation_weights(column_map)
    pcawg_tumor_tissues = get_pcawg_tumor_tissues()
    avg_t_per_chunk, chunk_size, n  = chunk2df(
        DATABASE,
        TABLE,
        TABLE_BURDEN,
        subpopulation_weights,
        column_map,
        pcawg_tumor_tissues
    )

if __name__ == "__main__":
    main()
