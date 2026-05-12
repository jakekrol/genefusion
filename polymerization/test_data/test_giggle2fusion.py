#!/usr/bin/env python3
import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *
from polymerization.giggle2fusion import *
from polymerization.polymerization import *
import time
import tempfile
import pandas as pd

OUTDIR_G2F='g2f_output'
OUTDIR_AGG='g2f_agg'
WORKERS=45
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
df_evidence = giggle2fusion(
    df_merged, df_shard, outdir=OUTDIR_G2F, logdir='g2f_log', max_workers=WORKERS, timeout=60*60*2, bgzip=True,
    sample_clean_func=lambda x: os.path.basename(x),
    path_bedfile=path_tmp_bed,
    gene_col_idx=3,
    evidence_right_gene_col=3,
    evidence_sample_col=14,
    outfile_prefix='',
    burden=True,
    bedtools_bin='/data/jake/bedtools.static.binary',
    verbose=True
)

# df_evidence = pd.read_csv('g2f_log/giggle_evidence_output.tsv', sep='\t')
agg_evidence_by_category(OUTDIR_G2F, OUTDIR_AGG, df_shard)