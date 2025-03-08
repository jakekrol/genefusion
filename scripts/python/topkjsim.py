#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os,sys
import plotly.express as px

#pcawg
# d='/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/2025_03_06-sample_top_fusion'
#1kg
d='/data/jake/genefusion/data/2024_11_01-fusions-1kg/2025_03_06-sample_top_fusion'
f = glob.glob(d+'/*')
l=[]
# loop over all samples
print('sets')
for f_i in f:
    A = set()
    df=pd.read_csv(f_i,sep='\t')
    s = os.path.basename(f_i).split('.')[0]
    # get top 100 fusion set
    for i,row in df.iterrows():
        v1 =row['v1']
        v2 =row['v2']
        fus=str(v1)+'-'+str(v2)
        A.add(fus)
    l.append((s,A))
l = sorted(l,key=lambda x:x[0])
print(l)
# do jaccard sim matrix
n = len(f)
J = np.zeros((n,n))
print('sim')
for i in range(n):
    for j in range(n):
        A = l[i][1]
        B = l[j][1]
        J[i,j] = len(A.intersection(B))/len(A.union(B))
pd.DataFrame(J).to_csv('jaccard_sim_1kg.tsv',index=False,sep='\t')
print('plot')
fig = px.imshow(J, 
                labels=dict(x="Sample", y="Sample", color="Jaccard Similarity"),
                x=[f'Sample {s}' for s,_ in l],
                y=[f'Sample {s}' for s,_ in l],
                color_continuous_scale='viridis')
fig.update_layout(title='Jaccard Similarity Heatmap')
fig.write_html('jaccard_sim_heatmap_1kg.html')
