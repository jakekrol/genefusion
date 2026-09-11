#!/usr/bin/env python3

import pandas as pd

def cleanup_bed(bed):
    df = pd.read_csv(bed,sep="\t",header=None)
    df.columns = ['chromosome', 'start', 'end', 'gene_name', 'strand']
    valid_chromosomes = list(range(1,23)) + ['X', 'Y']
    valid_chromosomes = [str(c) for c in valid_chromosomes]
    mask = df['chromosome'].isin(valid_chromosomes)
    df = df[mask]
    df.to_csv(bed,sep="\t", header=False,index=False)

for bed in ['grch37.genes.bed', 'grch37.genes.sort.bed']:
    cleanup_bed(bed)

