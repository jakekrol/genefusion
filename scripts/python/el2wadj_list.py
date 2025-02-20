#!/usr/bin/env python3

import sys
import os
import time
import pandas as pd

# took about 0.16 sec for a 100k edge list 

n = 25374 # number of genes
args = sys.argv[1:]
print(args)
# el = args[0]
el = "NA19443.el"
# out = args[1]
out = 'test.adj'

t_0 = time.time()
adj_list = {}
# add gene nodes
for i in range(n+1):
    adj_list[i] = {}
with open(el, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        i = int(line[0])
        j = int(line[1])
        # update i's adj list
        if j in adj_list[i]:
            adj_list[i][j] += 1
        else:
            adj_list[i][j] = 1
        # update j's adj list
        if i in adj_list[j]:
            adj_list[j][i] += 1
        else:
            adj_list[j][i] = 1
# write adjacency list to file
with open(out, 'w') as f:
    for k, v in adj_list.items():
        f.write(str(k) + '\t' + str(v) + '\n')
print('Time elapsed: ', time.time() - t_0)

