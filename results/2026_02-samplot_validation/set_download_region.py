#!/usr/bin/env python3
import argparse
import pandas as pd
import os,sys

parser=argparse.ArgumentParser(description="Set download region for each sample")
parser.add_argument("-i","--input",help="Input file with sample information")
parser.add_argument("-o","--output",help="Output file with download region information")
parser.add_argument("--len_multiplier", type=float,default=2.0, help="Get -/+ len_multiplier * breakpoint_length on each side of breakpoint as download region")
args=parser.parse_args()

grch37_chr_lengths = {
    "1": 249250621,
    "2": 243199373,
    "3": 198022430,
    "4": 191154276,
    "5": 180915260,
    "6": 171115067,
    "7": 159138663,
    "8": 146364022,
    "9": 141213431,
    "10": 135534747,
    "11": 135006516,
    "12": 133851895,
    "13": 115169878,
    "14": 107349540,
    "15": 102531392,
    "16": 90354753,
    "17": 81195210,
    "18": 78077248,
    "19": 59128983,
    "20": 63025520,
    "21": 48129895,
    "22": 51304566,
    "X": 155270560,
    "Y": 59373566,
    "M": 16569
}


def main():
    df=pd.read_csv(args.input,sep="\t")
    df['sv_length']=df['right_breakpoint']-df['left_breakpoint']
    df['region_start']=df['left_breakpoint']-(args.len_multiplier*df['sv_length'])
    df['region_start']=df['region_start'].apply(lambda x: max(1, x)).astype(int)
    df['region_end']=df['right_breakpoint']+(args.len_multiplier*df['sv_length'])
    df['region_end']=df.apply(lambda row: min(row['region_end'], grch37_chr_lengths.get(str(row['right_gene_chrom']), row['region_end'])), axis=1).astype(int)
    df['region_string']=df['right_gene_chrom'].astype(str)+":"+df['region_start'].astype(str)+"-"+df['region_end'].astype(str)
    df.to_csv(args.output,sep="\t",index=False)

if __name__ == "__main__":
    main()