#!/usr/bin/env python3
import argparse
import gzip
import os
import pandas as pd
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument("--gtf", default="../../data/2026_06-grch37-transcripts/gencode.v19.annotation.gtf.gz")
parser.add_argument("-output",default="exon_intervals.tsv.gz")
args = parser.parse_args()

# columns, 0-indexed
COL_CHROM=0
COL_RECORD_TYPE=2
COL_START=3
COL_END=4
COL_INFO=8

def validate_args():
    assert os.path.exists(args.gtf)

def main():
    with gzip.open(args.gtf, "rt") as f:
        data = defaultdict(dict)
        for i,line in enumerate(f.readlines()):
            # if i % 100000 == 0:
            #     print(f"# line {i}")
            if line.startswith("#"):
                continue
            fields = line.rstrip().split('\t')
            record_type = fields[COL_RECORD_TYPE]
            if record_type != 'exon':
                continue
            chrom = fields[COL_CHROM]
            start = fields[COL_START]
            end = fields[COL_END]
            if end <= start:
                continue
            info = fields[COL_INFO]
            info_fields = info.split(';')
            gene_name = None
            transcript_name = None
            for i in info_fields:
                if 'gene_name' in i:
                    gene_name = i.strip().split(' ')[1].replace('"', '')
                if 'transcript_name' in i:
                    transcript_name = i.strip().split(' ')[1].replace('"', '')
            if not all([chrom, start, end, gene_name, transcript_name]):
                continue
            else:
                if not (transcript_name in data[gene_name].keys()):
                    data[gene_name][transcript_name] = [(start,end)]
                else:
                    data[gene_name][transcript_name].append( (start,end) )
                    # sort, gtf is already sorted, but here we are pedantic
                    data[gene_name][transcript_name] = sorted(data[gene_name][transcript_name], key=lambda x: x[0])
    df_data = []
    for gene, transcript_dict in data.items():
        for transcript_name, exon_intervals in transcript_dict.items():
            df_data.append(
                (
                    gene,
                    transcript_name,
                    exon_intervals
                )
            )
    df = pd.DataFrame(df_data, columns = ['gene', 'transcript_name', 'exon_intervals'])
    m=df.shape[0]
    df.drop_duplicates(subset=["gene", "transcript_name"],inplace=True)
    m_drop = df.shape[0]
    print(f"# dropped duplicates: {m-m_drop}")
    df.reset_index(drop=True,inplace=True)
    print("# na counts")
    print(df.isna().sum())
    df.to_csv(args.output,sep="\t", index=False)
        

if __name__ == "__main__":
    main()