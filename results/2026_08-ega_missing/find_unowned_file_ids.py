#!/usr/bin/env python3
import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--icgc_metadata", default='../../data/2024_08-icgc_legacy_locations/icgc25k-legacy-data-locations.no_index.tsv', help="Path to the ICGC metadata file")
parser.add_argument("--file_ids", default='file_ids.tsv', help="Path to the file containing file IDs")
parser.add_argument("--output", default='unowned_file_ids.tsv', help="Path to the output file for unowned file IDs")
args = parser.parse_args()

FORMAT='BAM'
EXPERIMENTAL_STRATEGIES = {'RNA-Seq', 'WGS'}
EGA_COLS = {'ega_dataset_id', 'ega_analysis_id', 'ega_file_id', 'ega_run_id'}
LOCATION='EGA'

def main():
	df_icgc = pd.read_csv(args.icgc_metadata, sep="\t")
	df_owned = pd.read_csv(args.file_ids, sep="\t")
	owned_file_ids = set(df_owned['file_id'])
	mask = (df_icgc['Format'] == FORMAT) & \
    	(df_icgc['Experimental_Strategy'].isin(EXPERIMENTAL_STRATEGIES) & \
        	(~df_icgc['File_ID'].isin(owned_file_ids))) & \
				(df_icgc['location'] == LOCATION)
	df_unowned = df_icgc[mask].reset_index(drop=True)
	df_unowned.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()