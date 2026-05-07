#!/usr/bin/env python3
import os,sys
from polymerization.io import *
from polymerization.stix2fusion import *
from polymerization.giggle2fusion import *
import time

# depends on giggle2fusion.test_merge_fusionset_bed2giggle and test_bedtools_intersect.py 

path_intersect='TMEM120A.giggle.swap.intersect.bed.gz'
outfile='TMEM120A.evidence.tsv'
intersect2evidence(path_intersect, outfile, right_gene_col=3, sample_column=14, bgzip=True)
