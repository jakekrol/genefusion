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
parser.add_argument("--score-table", default="score-neg_tissue_agnostic.tsv", help="Path to score table tsv file")
parser.add_argument("--outdir", default="roc_data-neg_tissue_agnostic")
parser.add_argument("--w_normals", default="0.0,0.25,0.5,0.75,1.0", help="Comma-separated list of w_normal values")
parser.add_argument("--thousg-only", action="store_true", help="Whether to also plot ROC curves for 1000G-only scoring (no normal evidence)")
parser.add_argument("--title", default="ROC tissue-agnostic negative scoring", help="Title for ROC plot")

W_NORMALS = [float(x) for x in parser.parse_args().w_normals.split(",")]


def build_roc_input(
    df: pd.DataFrame,
    pos_col: str,
    neg_col: str,
    outfile: Path,
) -> None:
    cols_keep = ["gene_left", "gene_right"]

    p = df[cols_keep + [pos_col]].copy()
    p.columns = cols_keep + ["score"]
    p["label"] = 1

    n = df[cols_keep + [neg_col]].copy()
    n.columns = cols_keep + ["score"]
    n["label"] = 0

    df_out = pd.concat([p, n], axis=0)
    df_out = df_out.sort_values("score", ascending=False)
    df_out.to_csv(outfile, sep="\t", index=False)

    return outfile

def main() -> None:
    args = parser.parse_args()

    assert os.path.isfile(args.roc_script), f"ROC script not found at {args.roc_script}"
    assert os.path.isfile(args.score_table), f"Score table not found at {args.score_table}"
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.score_table, sep="\t")

    roc_files: list[str] = []
    names: list[str] = []

    for w_normal in W_NORMALS:
        if not args.thousg_only:
            # Negative set 1 (off-tissue scoring)
            outfile = Path(f"{args.outdir}/roc_input-neg1-wnormal_{w_normal}.tsv")
            build_roc_input(
                df=df,
                pos_col=f"score_as_positive-wnormal_{w_normal}",
                neg_col=f"score_as_negative-wnormal_{w_normal}",
                outfile=outfile,
            )
            roc_files.append(str(outfile))
            names.append(f"Wnormal={w_normal}")


        # Negative set 1 + 1000g weighting
        outfile_w_1000g = Path(f"{args.outdir}/roc_input-neg1-wnormal_{w_normal}-w_1000g.tsv")
        build_roc_input(
            df=df,
            pos_col=f"score_as_positive-w_1000g-wnormal_{w_normal}",
            neg_col=f"score_as_negative-w_1000g-wnormal_{w_normal}",
            outfile=outfile_w_1000g,
        )
        roc_files.append(str(outfile_w_1000g))
        names.append(f"Wnormal={w_normal};1000G")


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
