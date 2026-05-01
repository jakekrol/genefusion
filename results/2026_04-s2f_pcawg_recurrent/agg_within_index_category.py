#!/usr/bin/env python

from polymerization.io import *
from polymerization.stix2fusion import *
import time
import argparse
import os,sys

parser = argparse.ArgumentParser()
parser.add_argument("--shardfile",
					default='shardfile.tsv'
)
parser.add_argument("--outdir",
                    help='the outdir of s2f',
                    default='stix_output'
)
parser.add_argument('--outdir_agg',
                    default='stix_output_agg',
					help='the outdir for aggregated evidence')

args = parser.parse_args()
df_shard=read_stix_shardfile(args.shardfile)
agg_stix_evidence_by_category(args.outdir, args.outdir_agg, df_shard)
