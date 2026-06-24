#!/usr/bin/env python
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Extract gene names, transcript ids, and exon ids from GTF file.")
parser.add_argument("--input", type=str, required=True, help="Input GTF file without header.")
parser.add_argument("--output", type=str, required=True, help="Output file name for the extracted data.")
parser.add_argument("--extract_gene_name", action="store_true", help="Extract gene names.")
parser.add_argument("--extract_transcript_id", action="store_true", help="Extract transcript ids.")
parser.add_argument("--extract_exon_id", action="store_true", help="Extract exon ids.")
args = parser.parse_args()


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

def main():
    df = pd.read_csv(args.input, sep="\t", header=None)
    if args.extract_gene_name:
        df["gene_name"] = df[4].apply(extract_gene_name)
    if args.extract_exon_id:
        df["exon_id"] = df[4].apply(extract_exon_id)
    if args.extract_transcript_id:
        df["transcript_id"] = df[4].apply(extract_transcript_id)
    # df = df[[0, 1, 2, "gene_name", "exon_id", "transcript_id", 3]]
    if args.extract_gene_name and args.extract_exon_id and args.extract_transcript_id:
        df = df[[0, 1, 2, "gene_name", "exon_id", "transcript_id", 3]]
    elif args.extract_gene_name and args.extract_exon_id:
        df = df[[0, 1, 2, "gene_name", "exon_id", 3]]
    elif args.extract_gene_name and args.extract_transcript_id:
        df = df[[0, 1, 2, "gene_name", "transcript_id", 3]]
    elif args.extract_exon_id and args.extract_transcript_id:
        df = df[[0, 1, 2, "exon_id", "transcript_id", 3]]
    elif args.extract_gene_name:
        df = df[[0, 1, 2, "gene_name", 3]]
    elif args.extract_exon_id:
        df = df[[0, 1, 2, "exon_id", 3]]
    elif args.extract_transcript_id:
        df = df[[0, 1, 2, "transcript_id", 3]]
    else:
        raise ValueError("At least one of --extract_gene_name, --extract_exon_id, or --extract_transcript_id must be specified.")

    df.to_csv(args.output, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()