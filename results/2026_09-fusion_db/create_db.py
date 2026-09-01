#!/usr/bin/env python
import subprocess
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('--db', default='fusion.db')
args = parser.parse_args()

datasets = [
	### thousg
    (
        "thousg_rna_mage",
        "../2026_07-1kg-rna/g2f_agg/mage_short_read_1000g_rna-fusion_evidence.tsv",
    ),
	(
		"thousg_dna_low_coverage",
		"../2026_06-g2f-all_gene_pairs/g2f_agg/low_coverage_1000g_dna-fusion_evidence.tsv"
	),
	(
		"thousg_dna_high_coverage",
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


result = subprocess.run(
    ["sqlite3", "--version"],
    capture_output=True,
    text=True,
)

print(result.stdout)
print(result.stderr)

for name, path in datasets:
    assert os.path.exists(path)
    build_table(args.db, name, path)