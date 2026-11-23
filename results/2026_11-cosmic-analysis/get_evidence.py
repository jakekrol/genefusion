#!/usr/bin/env python3
import argparse
import pandas as pd
import sqlite3
parser=argparse.ArgumentParser()
# parser.add_argument("--fusions", default="cosmic_fusions.subset.tsv")
parser.add_argument("--fusions", default="cosmic_fusions.tsv")
parser.add_argument("--db", default="../2026_09-fusion_db/fusion.db")
parser.add_argument("--output",default="cosmic_fusions.evidence.tsv")
args = parser.parse_args()

TABLE_EVIDENCE='fusion_evidence'

def main():
    con = sqlite3.connect(args.db)
    df_fusion = pd.read_csv(args.fusions,sep="\t")
    df_fusion.to_sql('tmp', con, if_exists='replace', index=False)
    result = pd.read_sql_query(f"""
        SELECT tmp.*, {TABLE_EVIDENCE}.*
        FROM tmp
        LEFT JOIN {TABLE_EVIDENCE}
        USING (gene_left, gene_right)
    """, con)
    result.to_csv(args.output,sep='\t',index=False)


if __name__=="__main__":
    main()