#!/usr/bin/env python

# took about 15m

from polymerization.io import *
from polymerization.stix2fusion import *
import time
import argparse
import os,sys

parser = argparse.ArgumentParser()
parser.add_argument("--shardfile",
                    default='shardfile.tsv',
                    help="shardfile")
parser.add_argument("--outdir",
                    default='stix_output',
                    help="Output directory")
parser.add_argument("--fusionset",
                    default="recurrent_normal_fusions.tsv",
                    help="Fusion set file")
parser.add_argument("--bedfile",
					default='grch37.genes.bed.added',
                    help="BED file with gene metadata")
parser.add_argument("--cpus",
					default=30,
					type=int,
					help="Number of CPUs to use for parallel processing")
args = parser.parse_args()

# inputs
df_fusionset=read_fusion_set(args.fusionset)
df_bed=read_bed(args.bedfile, gene_col_idx=3)
df_shard=read_stix_shardfile(args.shardfile)
# setup
df_fusionset_sort = left_sort_fusion_set(df_fusionset, df_bed)
df_merged = merge_fusion_set_with_bed(df_fusionset_sort, df_bed)
# run
t_0=time.time()
merge_fusion_set_bed2stix(
	df_merged,
	df_shard,
	args.outdir,
	max_workers=args.cpus,
)
t_1=time.time()
print(f"Total time: {t_1-t_0} seconds")

