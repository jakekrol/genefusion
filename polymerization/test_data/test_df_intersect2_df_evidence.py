#!/usr/bin/env python3
import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *
from polymerization.giggle2fusion import *
from polymerization.polymerization import *
import time
import tempfile

OUTDIR='g2f_output'
WORKERS=30
df_fusionset=read_fusion_set('normal_fusion.modif.tsv')
df_bed=read_bed('grch37.genes.bed.added', gene_col_idx=3)
df_bed_subset = subset_bed_by_fusion_set(df_fusionset, df_bed)
# write subset bed to tmpfile
with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_bed:
    df_bed_subset.to_csv(tmp_bed.name, sep='\t', index=False, header=False)
    path_tmp_bed = tmp_bed.name

df_fusionset_sort = left_sort_fusion_set(df_fusionset, df_bed_subset)
df_merged = merge_fusion_set_with_bed(df_fusionset_sort, df_bed_subset)
df_shard=read_giggle_shardfile('shardfile.giggle.2.tsv')
# giggle queries
df_giggle = merge_fusion_set_bed2giggle(
    df_merged, df_shard, outdir=OUTDIR, max_workers=WORKERS, timeout=60*60*2, bgzip=True
)
print(df_giggle.head())

# swap intervals in giggle output
df_swap = giggle2swap(df_giggle, df_shard, bgzip=True)
print(df_swap.head())

# intersect
df_intersect = swap2intersect(
    df_swap,
    df_shard,
    path_bedfile=path_tmp_bed, # use the subset bed file for intersection
    gene_col_idx=3,
    bgzip=True,
    bedtools_bin='/data/jake/bedtools.static.binary',
    max_workers=WORKERS
)
print(df_intersect.head())

# evidence
df_evidence = df_intersect2df_evidence(
    df_intersect,
    df_shard,
    right_gene_col=3,
    sample_column=14,
    bgzip=True,
    burden=True,
    max_workers=WORKERS,
)
print(df_evidence.head())

# aggregate evidence
outdir_agg='g2f_agg'
agg_evidence_by_category(
    outdir_g2f=OUTDIR,
    outdir_agg=outdir_agg,
    df_giggle_shards=df_shard,
    outfile_suffix='-fusion_evidence.tsv',
    evidence_pattern='*.evidence.tsv'
)
breakpoint()
