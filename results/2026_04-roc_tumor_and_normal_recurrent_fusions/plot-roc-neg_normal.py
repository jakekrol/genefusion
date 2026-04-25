#!/usr/bin/env python3

import subprocess
from pathlib import Path
import argparse
import pandas as pd
import os, sys
import matplotlib.pyplot as plt


parser=argparse.ArgumentParser()
parser.add_argument("--roc-script", default="/data/jake/rl-tools/plot/roc.py",
                    help="Path to roc plotting script")
parser.add_argument("--score-table", default="score-all-tumor-normal-recurrent.tsv", help="Path to score table tsv file")
parser.add_argument("--outdir", default="roc_data-neg_normal")
parser.add_argument("--w_normals", default="0.0,0.25,0.5,0.75,1.0", help="Comma-separated list of w_normal values")
parser.add_argument("--thousg-only", action="store_true", help="Whether to also plot ROC curves for 1000G-only scoring (no normal evidence)")
parser.add_argument("--title", default="", help="Title for ROC plot")

W_NORMALS = [float(x) for x in parser.parse_args().w_normals.split(",")]


def main() -> None:
    args = parser.parse_args()

    assert os.path.isfile(args.roc_script), f"ROC script not found at {args.roc_script}"
    assert os.path.isfile(args.score_table), f"Score table not found at {args.score_table}"
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.score_table, sep="\t")
    score_cols = [col for col in df.columns if col.startswith("score")]


    roc_files: list[str] = []
    names: list[str] = []
    for w_normal in W_NORMALS:
        for score_col in score_cols:
            if (args.thousg_only) and ("w_1000g" not in score_col):
                continue
            if f"wnormal_{w_normal}" in score_col:
                if "w_1000g" in score_col:
                    outfile = Path(args.outdir) / f"roc_input-neg_normal-w_{w_normal}-w_1000g.tsv"
                    if args.thousg_only:
                        name = f"Wnormal={w_normal}"
                    else:
                        name = f"Wnormal={w_normal};1000G"
                else:
                    outfile = Path(args.outdir) / f"roc_input-neg_normal-w_{w_normal}.tsv"
                    name = f"Wnormal={w_normal}"
                df_subset = df[["gene_left", "gene_right", score_col,"label"]]
                # drop na
                df_subset = df_subset.dropna(subset=[score_col])
                df_subset.to_csv(outfile, sep="\t", index=False)
                roc_files.append(str(outfile))
                names.append(name)

    print("ROC files:", ",".join(roc_files))

    cmd = [
        str(Path(args.roc_script).resolve()),
        "--scores",
        ",".join(roc_files),
        "--output",
        f"{args.outdir}/roc.png",
        "--header",
        "--score_col",
        "2",
        "--label_col",
        "3",
        "--names",
        ",".join(names),
        "--title",
        f"{args.title}",
    ]

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
