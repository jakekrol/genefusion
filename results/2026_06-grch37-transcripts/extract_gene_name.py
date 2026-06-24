#!/usr/bin/env python
import pandas as pd

df = pd.read_csv("gencode.v19.annotation.gtf.nohead.exon", sep="\t", header=None)

def extract_gene_name(info_str):
    fields = info_str.split(";")
    for field in fields:
        if field.strip().startswith("gene_name"):
            return field.strip().split('"')[1]
    return None

def extract_exon_id(info_str):
    fields = info_str.split(";")
    for field in fields:
        if field.strip().startswith("exon_id"):
            return field.strip().split('"')[1]
    return None

def extract_transcript_id(info_str):
    fields = info_str.split(";")
    for field in fields:
        if field.strip().startswith("transcript_id"):
            return field.strip().split('"')[1]
    return None

df["gene_name"] = df[4].apply(extract_gene_name)
df["exon_id"] = df[4].apply(extract_exon_id)
df["transcript_id"] = df[4].apply(extract_transcript_id)
df = df[[0, 1, 2, "gene_name", "exon_id", "transcript_id", 3]]

df.to_csv("gencode.v19.annotation.gtf.exons.bed", sep="\t", index=False, header=False)
