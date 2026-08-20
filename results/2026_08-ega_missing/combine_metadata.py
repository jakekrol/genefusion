#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--dir_metadata", default="ega_metadata")
parser.add_argument("--merged_metadata_suffix", default="merged_metadata.tsv")
parser.add_argument("--output", default="combined_ega_metadata.tsv.gz")
args = parser.parse_args()

def main():
	metadata_files = glob.glob(os.path.join(args.dir_metadata,"*", f"*{args.merged_metadata_suffix}"))
	print(f"# found {len(metadata_files)} metadata files:")
	df = pd.DataFrame()
	for f in metadata_files:
		if os.path.getsize(f) == 0:
			print(f"# skipping empty file: {f}")
			continue
		df = pd.concat([df, pd.read_csv(f, sep="\t")], ignore_index=True)
	combined_df = df.drop_duplicates().reset_index(drop=True)
	combined_df.to_csv(args.output, sep="\t", index=False, compression="gzip")

if __name__ == "__main__":
	main()