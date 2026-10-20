#!/usr/bin/env python3
import argparse
import pandas as pd
import sqlite3
parser=argparse.ArgumentParser()
parser.add_argument("--db", default="../2026_09-fusion_db/fusion.db")
parser.add_argument("--output",default="evidence.tsv")
args = parser.parse_args()

TBL_EVIDENCE="fusion_evidence"
GENE_LEFT="PARGP1"
GENE_RIGHT="TIMM23"

def main():
	conn = sqlite3.connect(args.db)
	query = f"""
	SELECT * FROM {TBL_EVIDENCE} 
	WHERE {TBL_EVIDENCE}.gene_left = '{GENE_LEFT}' AND {TBL_EVIDENCE}.gene_right = '{GENE_RIGHT}';
	"""
	df = pd.read_sql_query(query, conn)
	df = df.T
	df = df.reset_index()
	df.columns = ['column', 'value']
	df.to_csv(args.output,sep="\t",index=False)
 

if __name__ == "__main__":
    main()