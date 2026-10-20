#!/usr/bin/env python3
import argparse
import os
import pandas as pd
from polymerization.analysis import intersect_file2breakpoints

GENE_LEFT="RP11-345J4.5"
GENE_RIGHT="RP11-345J4.8"
COL_SAMPLE='sample'
K=3

parser = argparse.ArgumentParser()
parser.add_argument("--pcawg_kidney_normal",default=f"../2026_06-g2f-all_gene_pairs/g2f_out/kidney_normal_rna/{GENE_LEFT}.giggle.clean.swap.intersect.bed.gz")
parser.add_argument("--pcawg2ega", default="../2026_10-pcawg2ega/pcawg_fileid2ega_fileid.tsv")
parser.add_argument("--output",default="samples.tsv")
args = parser.parse_args()


def main():
	df_pcawg2ega = pd.read_csv(args.pcawg2ega,sep="\t")
	df_bp = intersect_file2breakpoints(args.pcawg_kidney_normal, GENE_RIGHT, False, True)
	df_bp['File_ID'] = df_bp['sample'].apply(lambda x: os.path.basename(x).split('.')[0])
	df_merge = pd.merge(df_bp,df_pcawg2ega,on='File_ID',how='left')
	ser = df_bp[COL_SAMPLE].value_counts().head(K).to_frame().reset_index()
	df_merge = pd.merge(ser,df_merge,on='sample',how='left')
	cols_keep=['sample', 'count', 'ega_file_id']
	df_merge = df_merge[cols_keep]
	df_merge.drop_duplicates(inplace=True)
	df_merge.to_csv(args.output,sep="\t",index=False)

if __name__ == "__main__":
    main()

