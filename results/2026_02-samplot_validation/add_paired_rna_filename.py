#!/usr/bin/env python3

import pandas as pd
import argparse
import os,sys

parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('--donor_id_map', type=str,
    default='../2025_12-pcawg_donor_fid_map/icgc25k-wgs_rna-file_id-file_name_donor_map.tsv')
args = parser.parse_args()

df_in = pd.read_csv(args.input, sep='\t')
donor_id_map = pd.read_csv(args.donor_id_map, sep='\t')
donor_id_map = donor_id_map[donor_id_map['data_type'] == 'RNA-Seq']
donor_id_map.rename(columns={'file_name': 'rna_file_name', 'file_id': 'rna_file_id'}, inplace=True)
donor_id_map.rename(columns={'specimen_type': 'specimen_type_rna'}, inplace=True)
donor_id_map = donor_id_map[donor_id_map['rna_file_name'].str.contains('.STAR')]
donor_id_map = donor_id_map[['donor_id', 'rna_file_id', 'rna_file_name', 'specimen_type_rna']]
joint = df_in.merge(donor_id_map, left_on='ICGC_Donor', right_on='donor_id', how='left')
# filter out any rows where we don't have a paired RNA file name
joint = joint[~joint['rna_file_name'].isna()]
# only write if the dataframe is non-empty
if joint.shape[0] > 0:
    joint.to_csv(args.output, sep='\t', index=False)
else:
    print(f'No rows with paired RNA file name found in {args.input}. No output written to {args.output}.', file=sys.stderr)