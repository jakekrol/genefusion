#!/usr/bin/env python3
import argparse
import glob
import os
import pandas as pd
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("--dir_index_pcawg_dna", type=str, default='/data/jake/stix-pcawg-dna', help="Directory containing the index of PCAWG DNA data")
parser.add_argument("--dir_index_pcawg_rna", type=str, default='/data/jake/stix-pcawg-rna', help="Directory containing the index of PCAWG RNA data")
parser.add_argument("--metadata_suffix", type=str, default=".yaml")
parser.add_argument("--output", type=str, default="file_ids.tsv", help="Output file to write the file IDs")
args = parser.parse_args()

def metadata2fileids(metadata_file):
	with open(metadata_file, 'r') as f:
		metadata = yaml.safe_load(f)
	file_ids = []
	for _, filenames in metadata.items():
		for filename in filenames:
			filename = os.path.basename(filename)
			file_id = filename.split('.')[0]
			file_ids.append(file_id)
	return file_ids

def main():
	metadata_dna_index = glob.glob(os.path.join(args.dir_index_pcawg_dna, f"*{args.metadata_suffix}"))
	metadata_rna_index = glob.glob(os.path.join(args.dir_index_pcawg_rna, f"*{args.metadata_suffix}"))
	dna_file_ids = []
	rna_file_ids = []
	for f in metadata_dna_index:
		print(f"# dna metadata file: {f}")
		dna_file_ids.extend(metadata2fileids(f))
	for f in metadata_rna_index:
		print(f"# rna metadata file: {f}")
		rna_file_ids.extend(metadata2fileids(f))
	print(f"# dna file ids: {len(dna_file_ids)}")
	print(f"# rna file ids: {len(rna_file_ids)}")
	dna_file_ids.sort()
	rna_file_ids.sort()
	with open(args.output, 'w') as f:
		f.write("file_id\ttechnology\n")
		for file_id in dna_file_ids:
			f.write(f"{file_id}\tdna\n")
		for file_id in rna_file_ids:
			f.write(f"{file_id}\trna\n")


if __name__ == "__main__":
    main()


