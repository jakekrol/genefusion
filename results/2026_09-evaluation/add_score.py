#!/usr/bin/env python3
import argparse
import pandas as pd
import sqlite3
parser = argparse.ArgumentParser()
parser.add_argument("--db",default="../2026_09-fusion_db/fusion.db")
parser.add_argument("--eval", default="./fusion_eval.tsv")
parser.add_argument("--output", default="./fusion_eval.scored.tsv")
parser.add_argument("--table_name",default="score")
args = parser.parse_args()

conn = sqlite3.connect(args.db)
df = pd.read_csv(args.eval,sep="\t")
df.to_sql("temp", conn, if_exists="replace", index=False)
sql_cols = pd.read_sql_query(f"PRAGMA table_info({args.table_name})", conn)["name"]
sql_cols = [f"{args.table_name}.{c}" for c in sql_cols if c not in ("gene_left", "gene_right")]

query = f"""
    SELECT temp.*, {", ".join(sql_cols)}
    FROM temp
    LEFT JOIN {args.table_name}
      ON temp.gene_left = {args.table_name}.gene_left
     AND temp.gene_right = {args.table_name}.gene_right
"""
result = pd.read_sql_query(query, conn)
result.to_csv(args.output,sep="\t")