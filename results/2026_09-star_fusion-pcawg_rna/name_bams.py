#!/usr/bin/env python3

import argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("--ega", default="ega_bams.tsv")
parser.add_argument("--pcawg", default="icgc25-legacy-data.bam.tsv")
parser.add_argument("--output", default="ega_bams.named.tsv")
args = parser.parse_args()

def main():
    df_ega = pd.read_csv(args.ega,sep="\t")
    df_ega['key'] = df_ega['file_name'].apply(lambda x: x.split('/')[2].replace(".cip", ""))
    df_pcawg = pd.read_csv(args.pcawg,sep="\t",header=None)
    df_pcawg = df_pcawg[[1,2]]
    df_pcawg.columns = ['pcawg_file_id', 'pcawg_file_name']
    pcawg_filename2id = dict()
    for i,row in df_pcawg.iterrows():
        pcawg_filename2id[row['pcawg_file_name']] = row['pcawg_file_id']
    df_ega['pcawg_file_id'] = df_ega['key'].map(pcawg_filename2id)
    df_ega.to_csv(args.output,sep="\t",index=False)

if __name__ == "__main__":
    main()
