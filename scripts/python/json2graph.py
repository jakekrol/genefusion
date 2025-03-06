#!/usr/bin/env python3
import sys,os
import json
import networkx as nx
import pandas as pd
import time

j = 'test.json'
index='/data/jake/genefusion/data/genes.index'
with open(j, 'r') as f:
    aj = json.load(f)
self_loops=False
s_idx = pd.read_csv(index, sep='\t', header=None, usecols=[1]).squeeze() # squeeze into series
g = nx.Graph()
t0=time.time()
# note i,j need to be strings when accessing aj
for i in aj.keys():
    g.add_node(int(i), label=s_idx[int(i)])
    for j in aj[i]['edges'].keys():
        w = int(aj[i]['edges'][j])
        # handle self loops
        if int(j) == int(i):
            if self_loops:
                g.add_edge(int(i), int(j), weight=w)
            else:    
                continue
        # check if edge already stored since adj list is undirected
        # and we store both directions
        if g.has_edge(int(i), int(j)):
            continue
        else:
            g.add_edge(int(i), int(j), weight=w) 
print('Time to load graph:', time.time()-t0)
print(g)