#!/usr/bin/env python3
from genefusion.genefusion import json2graph
import networkx as nx
import json
import os,sys
f='FI9742.json'
out='test.json'
print('read')
g=json2graph(f, self_loops=True)
print('write')
def graph2json(g, out):
    aj={}
    for n in g.nodes():
        aj[n]={'edges':{}}
    for i,j,d in g.edges(data=True):
        w = d['weight']
        if j in aj[i]['edges'].keys():
            pass
        else:
            aj[i]['edges'][j] = w
        # update j's adj list
        if i in aj[j]['edges'].keys():
            pass
        else:
            aj[j]['edges'][i] = w
    with open(out, 'w') as f:
        json.dump(aj, f,indent=4,sort_keys=True)
graph2json(g, out)

    