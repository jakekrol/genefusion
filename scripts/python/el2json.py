#!/usr/bin/env python3

import sys,os
import json
el='/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/FI22097.el'
out = 'test.json'

n=25375
aj = {}
for i in range(0,n):
    aj[i] = {'edges':{}}


with open(el, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        i = int(line[0])
        j = int(line[1])
        # update i's adj list
        if j in aj[i]['edges'].keys():
            aj[i]['edges'][j] += 1
        else:
            aj[i]['edges'][j] = 1
        # update j's adj list
        if i in aj[j]['edges'].keys():
            aj[j]['edges'][i] += 1
        else:
            aj[j]['edges'][i] = 1
with open(out, 'w') as f:
    json.dump(aj, f, indent=4, sort_keys=True)
