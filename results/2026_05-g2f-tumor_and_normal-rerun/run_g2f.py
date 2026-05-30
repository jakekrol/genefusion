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

parser = argparse.ArgumentParser(description='Get gene fusion evidence for normal fusions')
parser.add_argument('--shardfile', default='./shardfile.tsv', help='path to giggle shardfile')
parser.add_argument('--bedfile', default='./grch37.genes.bed.added', help='path to bed file')
parser.add_argument('--fusion_set', default='./query_fusions.tsv', help='path to fusion set file')
parser.add_argument('--outdir_g2f', default='./g2f_out', help='output directory')
parser.add_argument('--outdir_agg', default='./g2f_agg', help='aggregated output directory')
parser.add_argument('--logdir', default='./g2f_log', help='giggle2fusion log directory')
parser.add_argument('--max_workers', type=int, default=55, help='number of parallel workers for giggle2fusion')
args = parser.parse_args()

os.makedirs(args.outdir_g2f, exist_ok=True)
os.makedirs(args.outdir_agg, exist_ok=True)
os.makedirs(args.logdir, exist_ok=True)

### read inputs
print(f"# reading giggle shard file: {args.shardfile}")
df_shard = read_giggle_shardfile(args.shardfile)
print(f"# reading bed file: {args.bedfile}")
df_bed = read_bed(args.bedfile, gene_col_idx=3)

print(f"# reading fusion set: {args.fusion_set}")
t_0=time.time()
df_fusions = read_fusion_set(args.fusion_set, reader='pandas')
print(f"# time to read fusion set: {time.time() - t_0:.2f} seconds")

print(f"# sorting fusion set by genomic position")
df_fusions_sort = left_sort_fusion_set(df_fusions, df_bed)

### merge fusion set with bed file
print(f"# merging fusion set with bed file")
# set columns in fusion set
df_fusions_sort.columns = ['gene_left', 'gene_right']
t_0=time.time()
df_merge = merge_fusion_set_with_bed(df_fusions_sort, df_bed, merger='pandas')
print(f"# time to merge fusion set with bed file: {time.time() - t_0:.2f} seconds")


# free memory
df_fusions=None
df_fusions_sort=None
df_bed=None

steps=['giggle', 'clean', 'swap', 'intersect', 'evidence']
df_evidence = giggle2fusion(
    df_merge, df_shard, outdir=args.outdir_g2f, logdir=args.logdir, max_workers=args.max_workers, timeout=60*60*2, bgzip=True,
    sample_clean_func=lambda x: os.path.basename(x),
    path_bedfile=args.bedfile,
    gene_col_idx=3,
    evidence_right_gene_col=3,
    evidence_sample_col=14,
    outfile_prefix='',
    burden=True,
    bedtools_bin='/data/jake/bedtools.static.binary',
    verbose=True,
    steps=steps
)

# free memory
df_merge=None

agg_evidence_by_category(args.outdir_g2f, args.outdir_agg, df_shard)
