#!/usr/bin/env python3

import pandas as pd
import argparse
import os,sys

parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('--donor_id_map', type=str,
    default='../results/2025_12-pcawg_donor_fid_map/icgc25k-wgs_rna-file_id-file_name_donor_map.tsv')
args = parser.parse_args()

df_in = pd.read_csv(args.input, sep='\t')
donor_id_map = pd.read_csv(args.donor_id_map, sep='\t')
breakpoint()