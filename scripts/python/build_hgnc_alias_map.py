#!/usr/bin/env python3
import gzip
import re
import pandas as pd
from collections import defaultdict
import argparse
import yaml
import os,sys

parser = argparse.ArgumentParser(description="Build HGNC alias table")
parser.add_argument("--hgnc", type=str, default="../../data/2025_09-gene_bed/genenames.tsv", help="Path to HGNC gene names TSV file")
parser.add_argument("--gtf", type=str, default="../../data/2025_09-gene_bed/gencode.v19.annotation.gtf.gz", help="Path to GENCODE v19 GTF file")
parser.add_argument("--output", "-o", type=str, default="hgnc_alias.yaml", help="Output YAML file name")
args = parser.parse_args()

# --- Load HGNC aliases ---
hgnc = pd.read_csv(args.hgnc, sep="\t", dtype=str)
hgnc.columns = hgnc.columns.str.strip()  # clean headers
print(hgnc.head())

def f(x):
    return x.upper() if type(x) == str else str(x).upper()

# gather all genes
# HGNC
X = set()
for _, row in hgnc.iterrows():
    a = row["Approved symbol"]
    a = f(a)
    X.add(a)
    if pd.notna(row.get("Alias symbols")):
        for b in row["Alias symbols"].split(", "):
            b = f(b)
            X.add(b)
    if pd.notna(row.get("Previous symbols")):
        for c in row["Previous symbols"].split(", "):
            c = f(c)
            X.add(c)

# GENCODE v19
with gzip.open(args.gtf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if cols[2] != "gene":
            continue
        attrs = cols[8]
        gnm = re.search(r'gene_name "([^"]+)"', attrs).group(1)
        gnm.replace("'", "")
        X.add(f(gnm))

# build alias map
data = {i: set() for i in X}
for _, row in hgnc.iterrows():
    a = row["Approved symbol"]
    a = f(a)
    if pd.notna(row.get("Alias symbols")):
        aliases = row["Alias symbols"].split(", ")
        for alias in aliases:
            c_alias = aliases.copy()
            c_alias.remove(alias)
            alias = f(alias)
            # add to canonical
            data[a].add(alias)
            # add aliases to each other, plus canonical
            for b in c_alias:
                b = f(b)
                data[alias].add(b)
                data[alias].add(a)
    if pd.notna(row.get("Previous symbols")):
        previous = row["Previous symbols"].split(", ")
        for prev in previous:
            c_prev = previous.copy()
            c_prev.remove(prev)
            prev = f(prev)
            # add to canonical
            data[a].add(prev)
            # add aliases to each other, plus canonical
            for b in c_prev:
                b  = f(b)
                data[prev].add(b)
                data[prev].add(a)
for k,v in data.items():
    v = list(data[k])
    v.sort()
    data[k] = v
sorted_dict = dict(sorted(data.items()))

with open(args.output, "w") as f:
    yaml.dump(sorted_dict,f)

print(f"Wrote {args.output} with {len(data)} entries")