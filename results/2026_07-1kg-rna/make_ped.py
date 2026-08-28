#!/usr/bin/env python

import argparse
import glob
from jkbiolib.datasets.loaders import *
import os
import pandas

parser = argparse.ArgumentParser()
parser.add_argument("--metadata", default="../../data/2026_08-thousg-rna-metadata/thousg-mage-sample-metadata.tsv")
parser.add_argument("--dir_bed", default="thousg_rna_bed_sort")
args = parser.parse_args()

def main():
	assert os.path.exists(args.metadata)
	assert os.path.isdir(args.dir_bed)
	df = pd.read_csv(args.metadata, sep="\t")
	bed_files = glob.glob(os.path.join(args.dir_bed), "*.bed.gz")

if __name__ == "__main__":
    main()