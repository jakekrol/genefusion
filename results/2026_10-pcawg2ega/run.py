#!/usr/bin/env python3
import argparse
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument("--icgc_metadata", default="../../data/2024_08-icgc_legacy_locations/icgc25k-legacy-data-locations.tsv.gz")
parser.add_argument("--output",default="pcawg_fileid2ega_fileid.tsv")
args = parser.parse_args()

def main():
    df = pd.read_csv(args.icgc_metadata,sep="\t",index_col=0)
    df=df[['File_ID','ega_file_id']]
    df = df.dropna()
    df.to_csv(args.output,sep="\t",index=False)

if __name__ == "__main__":
    main()