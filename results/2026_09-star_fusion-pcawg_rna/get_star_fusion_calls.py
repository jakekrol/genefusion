#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--dir_prefix", default="FI") # 1 dir per pcawg sample
parser.add_argument("--dir_star_fusion",default="star-fusion-out")
parser.add_argument("--output",default="thousg_rna-star_fusion_calls.tsv")
parser.add_argument("--calls_filename", default="star-fusion.fusion_predictions.tsv")
args = parser.parse_args()

def get_fusion_calls(path_star_fusion_pred):
    with open(path_star_fusion_pred,"r") as f:
        # skip header
        next(f)
        fusions=set()
        for line in f.readlines():
            # first field is fusion calls
            fusion = line.split('\t')[0].strip()
            fusions.add(fusion)
        return fusions


def main():
    files = glob.glob(f"{args.dir_prefix}*/{args.dir_star_fusion}/{args.calls_filename}") # FI*/star-fusion-outdir/star-fusion.fusion_predictions.tsv
    fusion_data = defaultdict(list)
    for path in files:
        p = Path(path)
        pcawg_file_id = p.parts[0]
        fusions = list(get_fusion_calls(path))
        fusion_data[pcawg_file_id] = fusions
    fusion_tbl_data = []
    for sample, fusions in fusion_data.items():
        for fusion in fusions:
            fusion_tbl_data.append((sample, fusion))
    df = pd.DataFrame(fusion_tbl_data, columns=['sample', 'fusion'])
    df.to_csv(args.output,sep="\t", index=False)

if __name__ == "__main__":
    main()