#!/usr/bin/env python3

import argparse
import pandas as pd
import os,sys

parser=argparse.ArgumentParser(description="Combine all plotting data into single table")
parser.add_argument("-i","--input_dir",help="Input directory with plotting data files")
parser.add_argument("-o","--output",help="Output file with combined plotting data")
args=parser.parse_args()

def main():

    files = os.listdir(args.input_dir)
    dfs=[]
    for f in files:
        df = pd.read_csv(os.path.join(args.input_dir, f), sep="\t")
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df = combined_df.drop_duplicates()
    combined_df = combined_df.sort_values(by=['tissue', 'svtype', 'tumor_paired_plus_split_samplewise_evidence'], ascending=[True, True, False])
    leading_cols = [
        'tissue', 'svtype', 'left', 'right', 'rna_file_id', 'region_string',
        'tumor_paired_plus_split_samplewise_evidence', 'ICGC_Donor', 'sample',
        'region_string', 'sv_length', 
        'left_breakpoint', 'right_breakpoint', 'left_gene_chrom', 'left_gene_start', 'left_gene_end', 
        'right_gene_chrom', 'right_gene_start', 'right_gene_end'
    ]
    other_cols = [col for col in combined_df.columns if col not in leading_cols]
    combined_df = combined_df[leading_cols + other_cols]    
    combined_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()