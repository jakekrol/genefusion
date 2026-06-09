#!/usr/bin/env python
import pandas as pd

df = pd.read_csv("gencode.v19.annotation.gtf.nohead.exon", sep="\t", header=None)

def extract_gene_name(info_str):
    fields = info_str.split(";")
    for field in fields:
        if field.strip().startswith("gene_name"):
            return field.strip().split('"')[1]
    return None

df["gene_name"] = df[4].apply(extract_gene_name)
df = df[[0, 1, 2, "gene_name", 3]]

df.to_csv("gencode.v19.annotation.gtf.exons.bed", sep="\t", index=False, header=False)
