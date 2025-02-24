#!/usr/bin/env python3
import os,sys
import pandas as pd

idx_erg=7172
idx_tmprss2=23014
args = sys.argv[1:]
ewdist= args[0]
# ewdist='/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/2025_02_20-ewdist/FI44231.ewdist'
# ewdist='/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/2025_02_20-ewdist/FI5528.ewdist'
# ewdist='/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/2025_02_20-ewdist/FI622.ewdist'
out= args[1]

df = pd.read_csv(ewdist, sep="\t")
# try both directions
x = df[df['node1']==idx_erg]
x = x[x['node2']==idx_tmprss2]
if len(x)!=0:
    w = x['weight'].values[0]
    with open(out, 'w') as f:
        f.write(str(w))
    sys.exit(0)
x = df[df['node1']==idx_tmprss2]
x = x[x['node2']==idx_erg]
if len(x)!=0:
    w=x['weight'].values[0]
    with open(out, 'w') as f:
        f.write(str(w))
    sys.exit(0)
else: 
    w = 0
    print('no fusion found')
with open(out, 'w') as f:
    f.write(str(w))