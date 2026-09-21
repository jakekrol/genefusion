#!/usr/bin/env python
import argparse
import os
import pandas as pd
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument('--db', default='fusion.db')
args = parser.parse_args()

TBL_GENE_PAIRS="sorted_gene_pairs"
TBL_FUSION_EVIDENCE="fusion_evidence"

datasets = [
	### thousg
    (
        "thousg_mage_short_read_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/mage_short_read_1000g_rna-fusion_evidence.tsv",
    ),
	(
		"thousg_low_coverage_dna",
		"../2026_06-g2f-all_gene_pairs/g2f_agg/low_coverage_1000g_dna-fusion_evidence.tsv"
	),
	(
		"thousg_high_coverage_dna",
		"../2026_06-g2f-all_gene_pairs/g2f_agg/high_coverage_1000g_dna-fusion_evidence.tsv"
	),
	### pcawg
	# blood
    (
        "pcawg_blood_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/blood_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_blood_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/blood_tumor_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_blood_tumor_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/blood_tumor_rna-fusion_evidence.tsv"
    ),
	# bone
    (
        "pcawg_bone_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/bone_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_bone_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/bone_tumor_dna-fusion_evidence.tsv"
    ),
	# breast
    (
        "pcawg_breast_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/breast_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_breast_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/breast_tumor_dna-fusion_evidence.tsv"
    ),
	# esophagus
    (
        "pcawg_esophagus_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/esophagus_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_esophagus_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/esophagus_tumor_dna-fusion_evidence.tsv"
    ),
	# gallbladder
    (
        "pcawg_gallbladder_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/gallbladder_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_gallbladder_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/gallbladder_tumor_dna-fusion_evidence.tsv"
    ),
	# headneck
    (
        "pcawg_headneck_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/headneck_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_headneck_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/headneck_tumor_dna-fusion_evidence.tsv"
    ),
	# kidney
    (
        "pcawg_kidney_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/kidney_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_kidney_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/kidney_tumor_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_kidney_normal_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/kidney_normal_rna-fusion_evidence.tsv"
    ),
    (
        "pcawg_kidney_tumor_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/kidney_tumor_rna-fusion_evidence.tsv"
    ),
	# liver
    (
        "pcawg_liver_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/liver_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_liver_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/liver_tumor_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_liver_normal_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/liver_normal_rna-fusion_evidence.tsv"
    ),
    (
        "pcawg_liver_tumor_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/liver_tumor_rna-fusion_evidence.tsv"
    ),
	# ovary
    (
        "pcawg_ovary_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/ovary_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_ovary_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/ovary_tumor_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_ovary_tumor_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/ovary_tumor_rna-fusion_evidence.tsv"
    ),
    # pancreas
    (
        "pcawg_pancreas_tumor_rna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/pancreas_tumor_rna-fusion_evidence.tsv"
    ),
	# prostate
    (
        "pcawg_prostate_normal_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/prostate_normal_dna-fusion_evidence.tsv"
    ),
    (
        "pcawg_prostate_tumor_dna",
        "../2026_06-g2f-all_gene_pairs/g2f_agg/prostate_tumor_dna-fusion_evidence.tsv"
    )
]

def build_table(db, name, path):
    sql = f"""
DROP TABLE IF EXISTS {name};
.print \"# creating {name} table\"
CREATE TABLE {name} (
    gene_left TEXT,
    gene_right TEXT,
    reads_{name} INTEGER,
    samples_{name} INTEGER
);

.mode tabs
.import --skip 1 '{path}' {name}

.print \"# indexing {name} table\"
CREATE UNIQUE INDEX idx_{name}
ON {name}(gene_left, gene_right);

.mode box
.headers on

.print \"# head of {name} table\"
SELECT * FROM {name} LIMIT 5;


.print \"# checking for duplicate fusions in {name} table\"
SELECT gene_left, gene_right, COUNT(*) AS n
FROM {name}
GROUP BY gene_left, gene_right
HAVING COUNT(*) > 1;
"""

    subprocess.run(
        ["sqlite3", db],
        input=sql,
        text=True,
        check=True,
    )

def union_fusion_keys(table_names, db, key1="gene_left", key2="gene_right", tbl_gene_pairs=TBL_GENE_PAIRS, tbl_fusion_evidence=TBL_FUSION_EVIDENCE, fill_value="0"):
    # delete table if already exists
    cmd = f"DROP TABLE IF EXISTS {joined_table_name};"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(
        ["sqlite3", db],
        input=cmd,
        text=True,
        check=True,
    )
    print(result.stdout)
    print(result.stderr)
    # get union of fusion keys
    cmd = f"CREATE TABLE {joined_table_name} AS WITH keys AS ( "
    n = len(table_names)
    for i,table in enumerate(table_names):
        if i < n-1:
            cmd += f"SELECT {key1}, {key2} FROM {table} UNION "
        else:
            cmd += f"SELECT {key1}, {key2} FROM {table}"
    cmd += ") SELECT * from keys;"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(
        ["sqlite3", db],
        input=cmd,
        text=True,
        check=True,
    )
    print(result.stdout)
    print(result.stderr)
    # index 
    cmd = f"CREATE INDEX idx_{joined_table_name} ON {joined_table_name}({key1}, {key2});"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(
        ["sqlite3", db],
        input=cmd,
        text=True,
        check=True,
    )
    print(result.stdout)
    print(result.stderr)

    # gather fusion evidence by join read and sample column data onto the union gene pairs table
    cmd = f"""
    DROP TABLE IF EXISTS {tbl_fusion_evidence};

    CREATE TABLE {tbl_fusion_evidence} AS
    SELECT
        {tbl_gene_pairs}.gene_left,
        {tbl_gene_pairs}.gene_right
    """

    for table in table_names:
        cmd += f"""
            , COALESCE({table}.reads_{table}, {fill_value}) AS reads_{table}
            , COALESCE({table}.samples_{table}, {fill_value}) AS samples_{table}
        """

    cmd += f"""
    FROM {tbl_gene_pairs}
    """

    for table in table_names:
        cmd += f"""
        LEFT JOIN {table}
            ON {tbl_gene_pairs}.gene_left = {table}.gene_left
            AND {tbl_gene_pairs}.gene_right = {table}.gene_right
        """

    cmd += ";"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(
        ["sqlite3", db],
        input=cmd,
        text=True,
        check=True,
    )
    print(result.stdout)
    print(result.stderr)

def add_burden(db,name, path, col_data_type_map):
    str_column_definition = ""
    for col, data_type in col_data_type_map.items():
        str_column_definition += f"{col} {data_type},\n"
    str_column_definition = str_column_definition[:-2]
        
    cmd = f"""
DROP TABLE IF EXISTS {name};
CREATE TABLE {name} (
{str_column_definition}
);
.mode tabs
.import --skip 1 '{path}' {name}
CREATE INDEX idx_{name} on {name}(gene);
.mode box
.print \"# preview {name} table\"
SELECT * FROM {name} LIMIT 5;
"""
    print(f"# running cmd: {cmd}")
    result = subprocess.run(
        ["sqlite3", db],
        input=cmd,
        text=True,
        check=True,
    )
    print(result.stdout)
    print(result.stderr)

    


### add individual datasets
print("# adding individual evidence tables")
for name, path in datasets:
    assert os.path.exists(path)
    build_table(args.db, name, path)

### join datasets
print("# joining fusion evidence tables")
table_names=[]
for name, _ in datasets:
    table_names.append(name)
union_fusion_keys(table_names, args.db)

### add burden table
print("# adding burden evidence table")
path="../2026_09-burden/burden.tsv"
name="burden"
df = pd.read_csv(path,sep='\t')
col_data_type_map = {}
for col in df.columns.tolist():
    if col == 'gene':
        col_data_type_map[col] = "TEXT"
    else:
        col_data_type_map[col] = "INTEGER"
add_burden(args.db, name, path, col_data_type_map)



