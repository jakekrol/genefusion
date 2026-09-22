#!/usr/bin/env python3
import argparse
import pandas as pd
import subprocess
import re
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("--ega_credentials", default="credentials.json", help="json config credential file")
parser.add_argument("--datasets", default="datasets-rna-star.txt", help="text file of ega dataset ids")
parser.add_argument("--output", default="ega_bams.tsv")
parser.add_argument("--ega_executable", default="egafetch-linux-amd64")
args = parser.parse_args()

# egalist returns a table-like structure, but the header is large and spaces are used instead of tabs for delim
def egalist2df(path: str):
    data = []
    table_start = 4 # 0-indexed, line 5
    with open(path, "r") as f:
        for i,line in enumerate(f.readlines()):
            # table start
            if i >= 4:
                # reached end of table
                if line == "\n":
                    break
                else:
                    # split by two or more spaces only
                    parts = re.split(r" {2,}", line)
                    parts = [p.strip() for p in parts]
                    data.append(parts)
    df = pd.DataFrame(data, columns = ['file_id', 'size', 'check', 'checksum', 'file_name'])
    return df
                    
        
    

def get_ega_dataset_files(id_dataset: str, ega_executable: str, credentials: str, outfile: str):
    cmd = f"{ega_executable} list {id_dataset} --cf {credentials} > {outfile}"
    print("# running cmd: {}".format(cmd))
    subprocess.run(
        cmd,
        shell=True,
        check=True
    )


def main():
    datasets = set()
    with open(args.datasets, 'r') as f:
        for line in f.readlines():
            datasets.add(line.strip())
    df_out = pd.DataFrame()
    for dataset in datasets:
        with tempfile.NamedTemporaryFile() as f:
            outfile = f.name
            get_ega_dataset_files(dataset, args.ega_executable, args.ega_credentials, outfile)
            df = egalist2df(outfile)
            df['is_bam'] = df['file_name'].apply(lambda x: x.endswith("bam.cip"))
            mask = df['is_bam'] == True
            df = df[mask].drop(columns=['is_bam']).reset_index(drop=True)
            df['dataset'] = dataset
            df_out = pd.concat([df_out, df], ignore_index=True)
    df_out.to_csv(args.output,sep="\t", index=False)

if __name__ == "__main__":
    main()
