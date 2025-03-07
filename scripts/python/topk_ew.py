#!/usr/bin/env python3
from genefusion.genefusion import json2graph
import networkx as nx
import os,sys
import heapq
import pandas as pd
args = sys.argv[1:]
j = args[0]
print('reading graph')
g = json2graph(j)

# topk edge weights
def topk_ew(G,k=100):
    mh = []
    for u, v, edge_data in g.edges(data=True):
        w = edge_data['weight']
        # if heap is not full, push the weight
        if len(mh) < k:
            heapq.heappush(mh, (w,u,v))

        # otherwise compare the weight with the smallest weight in the heap (element 0)
        else:
            if w > mh[0][0]:
                heapq.heapreplace(mh, (w,u,v))
    topk= sorted(mh,reverse=True, key=lambda x: x[0])
    return topk

print('topk edge weights')
topk = topk_ew(g)
df = pd.DataFrame(topk,columns=['weight','v1','v2'])
df.to_csv('test2.csv',index=False)
