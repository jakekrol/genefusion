#!/usr/bin/env bash
FILE=grch37.genes.sort.bed
sort -k1,1 -k2,2n -k3,3n $FILE > x && mv x $FILE