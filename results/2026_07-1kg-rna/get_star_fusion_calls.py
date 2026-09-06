#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--dir_star_fusion", default="star_fusion-parallel")
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
    files = glob.glob(f"{args.dir_star_fusion}/*/{args.calls_filename}")
    fusion_data = defaultdict(list)
    for path in files:
        sample = os.path.basename(os.path.dirname(path))
        fusions = list(get_fusion_calls(path))
        fusion_data[sample] = fusions
        # for fusion in fusions:
        #     fusion_counts
        #     is_fusion_in_sample_keys = fusion in fusion_counts[sample].keys()
        #     if not is_fusion_in_sample_keys:
        #         fusion_counts[sample][fusion] = 1
        #     else:
        #         fusion_counts[sample][fusion] +=1
    fusion_tbl_data = []
    for sample, fusions in fusion_data.items():
        for fusion in fusions:
            fusion_tbl_data.append((sample, fusion))
    df = pd.DataFrame(fusion_tbl_data, columns=['sample', 'fusion'])
    df.to_csv(args.output,sep="\t", index=False)

if __name__ == "__main__":
    main()