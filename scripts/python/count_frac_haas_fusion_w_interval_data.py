#!/usr/bin/env python3
import yaml
import argparse
import os,sys

parser = argparse.ArgumentParser(description="Count fraction of Haas fusions with interval data")
parser.add_argument("--fusions", type=str, default="./haas_fusions.txt", help="Input Haas fusions")
parser.add_argument("--genes", type=str, default="./genes_bed.txt", help="Input txt file of bed genes")
parser.add_argument("--alias", type=str, default="merged_alias.yaml", help="Input merged alias YAML file" )
parser.add_argument("--output", "-o", type=str, default="haas_fusion_interval_fraction.txt", help="Output text file")
args = parser.parse_args()

# parse haas fusions
fusions = set()
with open(args.fusions, 'r') as f:
    for line in f:
        fusions.add(line.strip().upper())
print(f"Read {len(fusions)} unique fusions from {args.fusions}")

# parse bed genes
genes = set()
with open(args.genes, 'r') as f:
    for line in f:
        genes.add(line.strip().upper())

print(f"Read {len(genes)} unique genes from {args.genes}")

# load merged alias map
with open(args.alias, 'r') as f:
    alias_map = yaml.safe_load(f)
print(f"Loaded {len(alias_map)} keys from merged alias map")

def search(query, genes, map):
    # search alias map
    if query in map.keys():
        key = query
    else:
        return False
    # see if we can route to a gene in bed file
    aliases = map[key]
    return any([i in genes for i in aliases + [key]])
        
        
n = len(fusions)
successes = 0
with open(args.output, 'w') as out:
    for i,x in enumerate(fusions):
        a,b = x.split('--')
        result = False
        if ',' in a:
            l_a = a.split(',')
            res_a = any([search(i,genes,alias_map) for i in l_a])
        else:
            res_a = search(a,genes,alias_map)
        if ',' in b:
            l_b = b.split(',')
            res_b = any([search(i,genes,alias_map) for i in l_b])
        else:
            res_b = search(b,genes,alias_map)
        result = res_a and res_b
        successes += result
        print(f"{i+1}/{n}: {x} -> {result}", file=out)
    print(successes/n, file=out)