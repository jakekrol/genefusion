#!/usr/bin/env bash
INFILE=grch37.genes.bed
OUTFILE=grch37.genes.sort.bed
sort -k1,1 -k2,2n -k3,3n $INFILE > $OUTFILE