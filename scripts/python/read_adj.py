#!/usr/bin/env python3

import sys
import os
import networkx as nx
import time
import ast
import pandas as pd

# took about 0.32 sec for a adj list with 100k edges

input = 'test.adj'
index='/data/jake/genefusion/data/genes.index'
s_idx = pd.read_csv(index, sep='\t', header=None, usecols=[1]).squeeze() # squeeze into series
t_0 = time.time()
G = nx.Graph()
with open(input, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        i = int(line[0])
        G.add_node(i)
        # add node label, could speed this up with binary search if necessary
        G.nodes[i]['label'] = s_idx[i]
        # dictionary of weighted edges
        d = ast.literal_eval(line[1].strip())
        for k, v in d.items():
            print(k)
            # check if edge already stored (since adj list is undirected)
            # if edge does not exist then skip
            if G.has_edge(i, k):
                continue
            else:
                G.add_edge(i, k, weight=v)
# print number of nodes and edges
print("###")
print(G.number_of_nodes())
print(G.number_of_edges())
print(G.nodes[25374])
# print the edges for node 25374
print(G.edges(25374))
print('Time elapsed: ', time.time() - t_0)
