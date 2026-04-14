#!/usr/bin/env python
from polymerization.io import *
from polymerization.stix2fusion import *

df_shard = read_stix_shardfile('shardfile.tsv')
print(df_shard.head())
