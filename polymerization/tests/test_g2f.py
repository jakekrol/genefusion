#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
from polymerization.giggle2fusion import *
from polymerization.io import *
from polymerization.polymerization import *
import time

SHARDFILE="input/shardfile.tsv"
BED="input/sampled_20.bed"
DIR_OUT="output/g2f_out"
DIR_AGG="output/g2f_agg"
DIR_LOG="output/g2f_log"

MAX_WORKERS=4

for directory in [DIR_OUT, DIR_AGG, DIR_LOG]:
    os.makedirs(directory, exist_ok=True)

### read inputs
df_shard = read_giggle_shardfile(SHARDFILE)
df_bed = read_bed(BED, gene_col_idx=3)

steps=['giggle', 'clean', 'swap', 'intersect', 'evidence']
df_evidence = giggle2fusion(
    df_bed,
    df_shard,
    DIR_OUT,
    DIR_LOG,
    BED,
    max_workers=MAX_WORKERS,
    bgzip=True,
    verbose=True,
    steps=['giggle', 'clean', 'swap', 'intersect', 'evidence']
)

agg_evidence_by_category(DIR_OUT,DIR_AGG, df_shard)

BED="input/sampled_21.bed"
df_bed = read_bed(BED, gene_col_idx=3)

df_evidence = giggle2fusion(
    df_bed,
    df_shard,
    DIR_OUT,
    DIR_LOG,
    BED,
    genes_only = set(["ADAM29"]),
    max_workers=MAX_WORKERS,
    bgzip=True,
    verbose=True,
    steps=['giggle', 'clean', 'swap', 'intersect', 'evidence']
)

agg_evidence_by_category(DIR_OUT,DIR_AGG, df_shard)

