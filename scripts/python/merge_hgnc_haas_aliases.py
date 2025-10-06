#!/usr/bin/env python3
import yaml
import argparse
import os,sys

parser = argparse.ArgumentParser(description="Merge HGNC and Haas alias tables")
parser.add_argument("--hgnc", type=str, default="hgnc_aliases.yaml", help="Input HGNC alias YAML file")
parser.add_argument("--haas", type=str, default="haas_aliases.yaml", help="Input Haas alias YAML file")
parser.add_argument("--output", type=str, default="merged_aliases.yaml", help="Output merged YAML file")
args = parser.parse_args()

print(f"Reading {args.hgnc} and {args.haas}...")
with open(args.hgnc, 'r') as f1, open(args.haas, 'r') as f2:
    hgnc_data = yaml.safe_load(f1)
    haas_data = yaml.safe_load(f2)

print(f"Merging aliases into {args.output}...")
merged = {}
for key in set(hgnc_data) | set(haas_data): # union of keys
    s1 = set(hgnc_data.get(key, []))
    s2 = set(haas_data.get(key, []))
    merged[key] = list(s1 | s2)

# Dump merged dictionary to YAML
print(f"Writing merged aliases to {args.output}...")
with open(args.output, "w") as f:
    yaml.dump(merged, f, sort_keys=True)