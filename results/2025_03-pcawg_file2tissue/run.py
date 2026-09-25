#!/usr/bin/env python3

import argparse
import pandas as pd
from polymerization.datasets import get_icgc_legacy_metadata

parser = argparse.ArgumentParser()
parser.add_argument("--project2tissue", default="../2024_09-icgc_project_abbreviations/icgc-abbrevs-simple.csv")
parser.add_argument("--output", default="pcawg_file_id2tissue.tsv")
args = parser.parse_args()

def main():
    df_proj2tissue = pd.read_csv(args.project2tissue)
    df_proj2tissue.columns = ['project', 'tissue']
    proj2tissue = dict()
    for row in df_proj2tissue.itertuples():
        proj2tissue[row[1]] = row[2]
    df_icgc_meta = get_icgc_legacy_metadata()
    df_icgc_meta['Project'] = df_icgc_meta['Project'].apply(lambda x: x.split('-')[0])
    df_icgc_meta['tissue'] = df_icgc_meta['Project'].map(proj2tissue)
    df_icgc_meta = df_icgc_meta[['File_ID', 'tissue']]
    df_icgc_meta.dropna(subset=['File_ID'],inplace=True)
    df_icgc_meta.fillna("unknown",inplace=True)
    df_icgc_meta.to_csv(args.output,sep="\t", index=False)

if __name__ == "__main__":
    main()