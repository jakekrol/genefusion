#!/usr/bin/env python

import argparse
import pandas as pd
parser=argparse.ArgumentParser()
parser.add_argument("--icgc_metadata", default="../../data/2024_08-icgc_legacy_locations/icgc25k-legacy-data-locations.no_index.tsv.gz")
parser.add_argument("--output", default="pcawg-ega-metadata-map.tsv.gz")
parser.add_argument("--project2tissue", default="../2024_09-icgc_project_abbreviations/icgc-abbrevs-simple.csv")
args = parser.parse_args()

COL_PCAWG_FILEID='File_ID'
COL_PCAWG_FILENAME='File_Name'
COL_PCAWG_DONOR='ICGC_Donor'
COL_PCAWG_SPECIMEN='Specimen_Type'
COL_PCAWG_FILETYPE='Data_Type'
COL_PCAWG_SEQTYPE='Experimental_Strategy'
COL_EGA_FILEID='ega_file_id'
COL_PCAWG_PROJECT='Project'
COL_PCAWG_SAMPLEID='ICGC_Sample_ID'
FILLNA='.'

def simple_specimen(x: str):
    x=x.lower()
    if ('tumour' in x) or ('tumor' in x):
        return 'tumor'
    if ('normal' in x):
        return 'normal'
    return x

def main():
    df_proj2tissue=pd.read_csv(args.project2tissue)
    df_proj2tissue.columns=[COL_PCAWG_PROJECT, 'tissue']
    s_proj2tissue=df_proj2tissue.set_index(COL_PCAWG_PROJECT).squeeze()
    df = pd.read_csv(args.icgc_metadata,sep="\t")
    df = df[
        [
            COL_PCAWG_FILEID,
            COL_PCAWG_DONOR,
            COL_PCAWG_SPECIMEN,
            COL_PCAWG_SAMPLEID,
            COL_PCAWG_FILETYPE,
            COL_PCAWG_SEQTYPE,
            COL_PCAWG_PROJECT,
            COL_EGA_FILEID
        ]
    ]
    df.fillna(FILLNA,inplace=True)
    mask=df[COL_PCAWG_FILETYPE].isin(['BAM', 'BAI'])
    df=df[mask]
    mask=df[COL_PCAWG_SEQTYPE].isin(['WGS','RNA-Seq'])
    df[COL_PCAWG_SPECIMEN] = df[COL_PCAWG_SPECIMEN].apply(simple_specimen)
    df[COL_PCAWG_PROJECT] = df[COL_PCAWG_PROJECT].apply(lambda x: x.split('-')[0])
    df['tissue'] = df[COL_PCAWG_PROJECT].map(s_proj2tissue)
    df=df.sort_values(COL_PCAWG_SAMPLEID)
    df.to_csv(args.output,index=False,sep="\t")

if __name__ == "__main__":
    main()

