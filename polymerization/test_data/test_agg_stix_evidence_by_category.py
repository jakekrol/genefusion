#!/usr/bin/env python3
import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *


df_shard=read_stix_shardfile('shardfile.tsv')

outdir='stix_output'
outdir_agg='stix_evidence_agg'
agg_stix_evidence_by_category(outdir, outdir_agg, df_shard)