#!/usr/bin/env python3
import argparse
import pandas as pd
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument("--unowned_icgc", default="unowned_file_ids.tsv")
parser.add_argument("--ega_metadata", default="combined_ega_metadata.tsv")
parser.add_argument("--output", default="unowned_icgc_and_ega_metadata.tsv")
args = parser.parse_args()

def main():
    df_icgc = pd.read_csv(args.unowned_icgc, sep="\t")
    df_ega = pd.read_csv(args.ega_metadata, sep="\t")
    # build a map of ega file ids to run ids
    ega_run_to_file = defaultdict(list)
    for i, row in df_ega.iterrows():
        if i % 10000 == 0:
            print(f"Processed {i} rows of EGA metadata")
        ega_run_to_file[row["run_accession_id"]].append(row["file_accession_id"])
    values_list_size = []
    for values_list in ega_run_to_file.values():
        values_list_size.append(len(values_list))
    breakpoint()
    # add a column to df_icgc for the corresponding ega
    # try matching file ids
    df_merge_fid = pd.merge(df_icgc, df_ega, how="inner", left_on="ega_file_id", right_on="file_accession_id")
    # try matching run ids
    df_merge_rid = pd.merge(df_icgc, df_ega, how="inner", left_on="ega_run_id", right_on="run_accession_id")
    breakpoint()
    df_merged.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()