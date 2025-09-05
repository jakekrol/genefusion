#!/usr/bin/env python3
import gzip
import re
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="Build a BED file with gene coordinates and aliases from HGNC and GENCODE.")
parser.add_argument("--hgnc", type=str, default="../../data/2025_09-gene_bed/genenames.tsv", help="Path to HGNC gene names TSV file")
parser.add_argument("--gtf", type=str, default="../../data/2025_09-gene_bed/gencode.v19.annotation.gtf.gz", help="Path to GENCODE v19 GTF file")
parser.add_argument("--output", "-o", type=str, default="grch37.chr.bed", help="Output BED file name")
args = parser.parse_args()

# --- Load HGNC aliases ---
hgnc = pd.read_csv(args.hgnc, sep="\t", dtype=str)
hgnc.columns = hgnc.columns.str.strip()  # clean headers
print(hgnc.head())

print("Building alias map...")
# Build alias → canonical map
alias_map = {}
for _, row in hgnc.iterrows():
    canon = str(row["Approved symbol"])
    alias_map[canon.upper()] = canon
    if pd.notna(row.get("Alias symbols")):
        for a in row["Alias symbols"].split(", "):
            alias_map[a.upper()] = canon
    if pd.notna(row.get("Previous symbols")):
        for p in row["Previous symbols"].split(", "):
            alias_map[p.upper()] = canon

# Precompute reverse mapping: canonical → list of aliases
reverse_map = defaultdict(list)
for alias, canon in alias_map.items():
    reverse_map[canon].append(alias)

# --- Load GENCODE v19 genes ---
print("Loading GENCODE v19 gene coordinates...")
genes = []
with gzip.open(args.gtf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if cols[2] != "gene":
            continue
        attrs = cols[8]
        gid = re.search(r'gene_id "([^"]+)"', attrs).group(1)
        gnm = re.search(r'gene_name "([^"]+)"', attrs).group(1)
        genes.append({
            "chr": cols[0],
            "start": int(cols[3]),
            "end": int(cols[4]),
            "strand": cols[6],
            "ensembl_gene_id": gid,
            "gene_name": gnm
        })
genes = pd.DataFrame(genes)

# ---------- Write to BED file ----------
print("Writing BED file with strand info...")
with open(args.output, "w") as out:
    for i, row in genes.iterrows():
        canon = row["gene_name"]
        chrom = row["chr"]
        start = row["start"] - 1  # BED is 0-based
        end = row["end"]
        strand = row["strand"]

        # Canonical name
        out.write(f"{chrom}\t{start}\t{end}\t{canon}\t.\t{strand}\n")

        # Aliases
        for alias in reverse_map.get(canon, []):
            if alias != canon.upper():
                out.write(f"{chrom}\t{start}\t{end}\t{alias}\t.\t{strand}\n")

        if (i + 1) % 1000 == 0 or i == len(genes) - 1:
            print(f"Processed {i+1}/{len(genes)} genes", end="\r")

print(f"\nBED file written to: {args.output}")
