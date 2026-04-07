#!/usr/bin/env python3
import os,sys
import pandas as pd
import argparse
import shutil
import time

t_0=time.time()
parser = argparse.ArgumentParser(description='split excord files into tumor/normal by filename')
parser.add_argument('-i', '--indir', type=str, required=True, help='Dir with excord files')
parser.add_argument('-l', '--lookup', type=str, default='/data/jake/genefusion/results/2025_05-pcawg_fileid2sample_type/fileid2sampletype.tsv', help='Path to file id to sample type mapping')
parser.add_argument('--outdir_normal', type=str, required=True, help='Normal outdir')
parser.add_argument('--outdir_tumor', type=str, required=True, help='Tumor outdir')
args = parser.parse_args()


def validation(args):
    assert os.path.isdir(args.indir), f'Input directory {args.indir} does not exist'
    assert os.path.isfile(args.lookup), f'Lookup file {args.lookup} does not exist'
    os.makedirs(args.outdir_normal, exist_ok=True)
    os.makedirs(args.outdir_tumor, exist_ok=True)

def main():
    print("# validating args")
    validation(args)
    files = [f for f in os.listdir(args.indir) if f.endswith('.excord.bed.gz')]
    files = [os.path.join(args.indir,f) for f in files]
    df_l = pd.read_csv(args.lookup, sep='\t')
    n = len(files) + 1
    for f in files:
        base = os.path.basename(f)
        fid = base.split('.')[0]
        specimen_type = df_l.loc[df_l['File_ID'] == fid, 'Specimen_Type'].values
        if len(specimen_type) == 0:
            raise ValueError(f'Could not find specimen type for file id {fid} from file {f}')
        specimen_type = specimen_type[0]
        if specimen_type.lower() == 'normal':
            shutil.copy(f, os.path.join(args.outdir_normal, base))
        elif specimen_type.lower() == 'tumour':
            shutil.copy(f, os.path.join(args.outdir_tumor, base))
        else:
            raise ValueError(f'Unknown specimen type {specimen_type} for file id {fid} from file {f}')
        n -= 1
        print(f'# copied {f}. {n-1} files remaining.')
    print(f"# done in {time.time()-t_0:.2f} seconds")

if __name__ == "__main__":
    main()
    
