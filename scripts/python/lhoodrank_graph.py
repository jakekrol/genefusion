#!/usr/bin/env python3
from genefusion.genefusion import json2graph
from genefusion.genefusion import graph2json
import pandas as pd
import numpy as np
from scipy.stats import nbinom
import os,sys
import math
args = sys.argv[1:]
j = args[0]
out = args[1]
# likelihood function
mu=120
theta=1.5
def nb_lhood(x, mu, theta):
    p = theta / (mu + theta)  # use success probability parameterization
    try: 
        lhood = math.log(nbinom.pmf(x, theta, p))
    except ValueError:
        lhood = -1000000
    return lhood
# read graph
g = json2graph(j)
# set weights to likelihoods
for i,j,d in g.edges(data=True):
    x = int(d['weight'])
    l = nb_lhood(x, mu, theta)
    d['weight'] = l # operates in place
# write graph
graph2json(g, out)
