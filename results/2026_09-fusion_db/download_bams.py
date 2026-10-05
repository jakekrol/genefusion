#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
import subprocess
parser=argparse.ArgumentParser()
parser.add_argument("--indir", default="top_fusions")
parser.add_argument("--outdir", default="bams")
parser.add_argument("--infile_suffix", default=".annotate.repeat_filt.breakpoint.rank.samples.tsv")
parser.add_argument("--egafetch_executable", default="egafetch-linux-amd64")
parser.add_argument("--ega_credentials",default="credentials.json")
parser.add_argument("--pcawg_file_id2ega_file_id",default="../2026_10-pcawg2ega/pcawg_fileid2ega_fileid.tsv")
parser.add_argument("--log", default="download_bams.log")
args = parser.parse_args()
os.makedirs(args.outdir,exist_ok=True)

def ega_download(file_id: str, path_credentials_json: str, executable: str, outdir: str):
    cmd = f"{executable} download {file_id} --output {outdir} --config-file {path_credentials_json} --restart"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(cmd, shell=True)
    return result.returncode

def main():
    df_pcawg2ega = pd.read_csv(args.pcawg_file_id2ega_file_id,sep="\t")
    df_pcawg2ega.set_index("File_ID",inplace=True)
    files = glob.glob(f"{args.indir}/*{args.infile_suffix}")
    ega_file_ids = set()
    # gather ids
    for f in files:
        df = pd.read_csv(f,sep="\t")
        top_samples = df['top_samples'].tolist()
        # parse pcawg file ids from csv strings
        for csv_string in top_samples:
            pcawg_file_ids = csv_string.split(",")
            for pcawg_file_id in pcawg_file_ids:
                try:
                    ega_file_id=df_pcawg2ega.loc[pcawg_file_id, "ega_file_id"]
                    ega_file_ids.add(ega_file_id)
                except KeyError:
                    print(f"# failed to lookup egaID for pcawgID={pcawg_file_id}")
                    continue
    # download
    results = []
    for ega_file_id in ega_file_ids:
        result = ega_download(ega_file_id, args.ega_credentials, args.egafetch_executable, args.outdir)
        results.append((ega_file_id,result))
    # log
    with open(args.log, "w") as f:
        for ega_file_id, result in results:
            f.write(f"{ega_file_id}\t{result}\n")
            

if __name__=="__main__":
    main()