#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
from polymerization.giggle2fusion import *
from polymerization.stix2fusion import *
from polymerization.io import *
from polymerization.polymerization import *
import time

parser = argparse.ArgumentParser(description='Aggregate fusion evidence')
parser.add_argument('--shardfile', default='./shardfile.tsv', help='path to giggle shardfile')
parser.add_argument('--outdir_g2f', default='./g2f_out', help='output directory')
parser.add_argument('--outdir_agg', default='./g2f_agg', help='aggregated output directory')
args = parser.parse_args()

os.makedirs(args.outdir_agg, exist_ok=True)

print(f"# reading giggle shard file: {args.shardfile}")
df_shard = read_giggle_shardfile(args.shardfile)

agg_evidence_by_category(args.outdir_g2f, args.outdir_agg, df_shard)