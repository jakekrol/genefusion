#!/usr/bin/env python3
from genefusion.genefusion import topk_ew
from genefusion.genefusion import json2graph
import os,sys
import pandas as pd
args = sys.argv[1:]
j = args[0]
out = args[1]
k=100
g = json2graph(j)
df = topk_ew(g,k)
df.to_csv(out,index=False, sep='\t')
