#!/usr/bin/env python3

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--fileid2tissue", default="../2025_03-pcawg_file2tissue/pcawg_file_id2tissue.tsv")
parser.add_argument("--fileid2specimen",default="../2025_05-pcawg_fileid2sample_type/fileid2sampletype.tsv")
parser.add_argument("--beds", default= "bed_clean_sort")
args = parser.parse_args()

def main():
    df_fid2tis = pd.read_csv(args.fileid2tissue, sep="\t")
    df_fid2tis.set_index("File_ID", inplace=True)
    df_fid2specimen = pd.read_csv(args.fileid2specimen,sep="\t")
    df_fid2specimen['Specimen_Type'] = df_fid2specimen['Specimen_Type'].replace('tumour', 'tumor')
    df_fid2specimen.set_index("File_ID", inplace=True)
    s_fid2specimen = df_fid2specimen['Specimen_Type'].squeeze()
    s_fid2tis = df_fid2tis['tissue'].squeeze()
    bed_files = os.listdir(args.beds)
    for f in bed_files:
        pcawg_file_id = f.split('.')[0]
        tissue = s_fid2tis[pcawg_file_id]
        specimen = s_fid2specimen[pcawg_file_id]
        dir_tissue = f"pcawg_rna_{tissue}_{specimen}_bed"
        os.makedirs(dir_tissue,exist_ok=True)
        destination = os.path.join(dir_tissue, f)
        src = os.path.join(args.beds,f)
        os.rename(src=src, dst=destination)

if __name__ == "__main__":
    main()
