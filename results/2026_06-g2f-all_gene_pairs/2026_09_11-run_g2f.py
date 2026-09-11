#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
from polymerization.giggle2fusion import *
from polymerization.io import *
from polymerization.polymerization import *
import time

SHARDFILE="shardfile.tsv"
BED="../2026_09-living_bed/grch37.genes.sort.bed"
DIR_OUT="g2f_out"
DIR_AGG="2026_09-g2f_agg"
DIR_LOG="g2f_log"
# FILE_GENES_TO_RUN="2026_09_11-genes_to_run.txt"

MAX_WORKERS=30

for directory in [DIR_OUT, DIR_AGG, DIR_LOG]:
    os.makedirs(directory, exist_ok=True)

### read inputs
df_shard = read_giggle_shardfile(SHARDFILE)
df_bed = read_bed(BED, gene_col_idx=3)
# genes_only = set()
# with open(FILE_GENES_TO_RUN, "r") as f:
#     for gene in f.readlines():
#         genes_only.add(gene.strip())

steps=['giggle', 'clean', 'swap', 'intersect', 'evidence']
df_evidence = giggle2fusion(
    df_bed,
    df_shard,
    DIR_OUT,
    DIR_LOG,
    BED,
    # genes_only=genes_only,
    max_workers=MAX_WORKERS,
    bgzip=True,
    verbose=True,
    steps=['giggle', 'clean', 'swap', 'intersect', 'evidence']
)

agg_evidence_by_category(DIR_OUT,DIR_AGG, df_shard)