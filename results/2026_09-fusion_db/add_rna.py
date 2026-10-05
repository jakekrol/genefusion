#!/usr/bin/env python3

import subprocess
import os

DIR_AGG='/data/jake/genefusion/results/2026_06-g2f-all_gene_pairs/g2f_agg'
RNA_TABLES=[
    ("blood_tumor_rna-fusion_evidence.tsv", "pcawg_blood_tumor_rna"),
    ("kidney_normal_rna-fusion_evidence.tsv", "pcawg_kidney_normal_rna"),
    ("kidney_tumor_rna-fusion_evidence.tsv", "pcawg_kidney_tumor_rna"),
    ("liver_normal_rna-fusion_evidence.tsv", "pcawg_liver_normal_rna"),
    ("liver_tumor_rna-fusion_evidence.tsv", "pcawg_liver_tumor_rna"),
    ("ovary_tumor_rna-fusion_evidence.tsv", "pcawg_ovary_tumor_rna"),
    ("pancreas_tumor_rna-fusion_evidence.tsv", "pcawg_pancreas_tumor_rna")
]

# def write_table(path, name,db ="fusion.db"):
#     sql = f"""
# DROP TABLE IF EXISTS {name};
# CREATE TABLE {name} (
#     gene_left TEXT,
#     gene_right TEXT,
#     reads_{name} INTEGER,
#     samples_{name} INTEGER
# );
# .mode tabs
# .import --skip 1 '{path}' {name}
# CREATE UNIQUE INDEX idx_{name}
# ON {name}(gene_left, gene_right);
# """
#     subprocess.run(
#         ["sqlite3", "--echo", db],
#         input=sql,
#         text=True,
#         check=True
#     )

# for path,name in RNA_TABLES:
#     path = os.path.join(DIR_AGG,path)
#     write_table(path,name)

def join(A_name, B_name, read_column_prefix="reads",sample_column_prefix="samples", key1="gene_left",key2="gene_right", fill_value=0):
	sql = f"""
UPDATE {A_name}
SET
    {read_column_prefix}_{B_name} = COALESCE((
        SELECT {B_name}.{read_column_prefix}_{B_name}
        FROM {B_name}
        WHERE {B_name}.{key1} = {A_name}.{key1} AND {B_name}.{key2} = {A_name}.{key2}
    ), {fill_value}),
    {sample_column_prefix}_{B_name} = COALESCE((
        SELECT {B_name}.{sample_column_prefix}_{B_name}
        FROM {B_name}
        WHERE {B_name}.{key1} = {A_name}.{key1} AND {B_name}.{key2} = {A_name}.{key2}
    ), {fill_value});
"""
	subprocess.run(
		["sqlite3", "--echo", "fusion.db"],
		input=sql,
		text=True,
		check=True
	)


A_name = "fusion_evidence"
B_names = [x[1] for x in RNA_TABLES]
for B_name in B_names:
    join(A_name,B_name)