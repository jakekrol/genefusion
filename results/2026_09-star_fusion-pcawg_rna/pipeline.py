#!/usr/bin/env python3

import argparse
import glob
import os
from pathlib import Path
import pandas as pd
import shutil
import subprocess
import time
from multiprocessing import Pool

parser=argparse.ArgumentParser()
parser.add_argument("--conda_env", default="star_fusion")
parser.add_argument("--ega_table", default="ega_bams.named.tsv")
parser.add_argument("--cpus_parallel", default=5, help="num. of input BAMs to process in parallel. keep low to minimize disk space usage")
parser.add_argument("--cpus_per_run", default=4, help="num of sub processes to be called from command line tools such as samtools and star fusion")
parser.add_argument("--ega_executable", default="egafetch-linux-amd64",help="https://github.com/khan-lab/EGAfetch")
parser.add_argument("--ega_credentials", default="credentials.json",help="path to ega credentials json")
parser.add_argument("--cache", action="store_true")
parser.add_argument("--genome_lib_dir", default="/data/jake/FusionAnnotator/genome_lib_dir", help="required by star fusion")
parser.add_argument("--outdir_star_fusion", default="star-fusion-out")
parser.add_argument("--logfile", default="pipeline.log")
parser.add_argument("--failed_rows", default="pipeline-failures.tsv", help="tsv of failed runs")
args = parser.parse_args()

def ega_download(file_id: str, path_credentials_json: str, executable: str, outdir: str):
    cmd = f"{executable} download {file_id} --output {outdir} --config-file {path_credentials_json}"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(cmd, shell=True)
    return result.returncode

def samtools_sort(bam_in: str, bam_out: str, threads: int = 4):
    cmd = f"samtools sort -@ {threads} -n -o {bam_out} {bam_in}"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(cmd, shell=True)
    return result.returncode

def bamtofastq(bam_sort: str, r1_out: str, r2_out: str):
    '''
    bam must be sorted by query name: samtools sort -n
    '''
    cmd = f"bedtools bamtofastq -i {bam_sort} -fq {r1_out} -fq2 {r2_out}"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(cmd, shell=True)
    return result.returncode

def star_fusion(fq_1: str, fq_2: str, outdir: str, star_fusion_conda: str, cpus: int, genome_lib_dir: str):
    cmd = "source $HOME/.bashrc; " # for conda init
    cmd += f"conda activate {star_fusion_conda}; "
    cmd += f"STAR-Fusion --CPU {cpus} --genome_lib_dir {genome_lib_dir} --left_fq {fq_1} --right_fq {fq_2} --output_dir {outdir}"
    print(f"# running cmd: {cmd}")
    result = subprocess.run(cmd, shell=True)
    return result.returncode

def validate_args():
    assert os.path.exists(args.ega_table)
    assert os.path.exists(args.ega_credentials)
    assert os.path.isdir(args.genome_lib_dir)
    assert shutil.which(args.ega_executable)

def validate_executables():
    assert shutil.which("samtools")
    assert shutil.which("bedtools")

def process_row(row):
    ega_file_id = row['file_id']
    pcawg_file_id = row['pcawg_file_id']
    os.makedirs(pcawg_file_id,exist_ok=True)
    ### download
    # try searching for cached bam
    bam_search =  glob.glob(f"{pcawg_file_id}/{ega_file_id}/*.bam")
    bam = bam_search[0] if bam_search else None
    if args.cache and bam:
        result_download = 0
    else:
        result_download = ega_download(
            file_id=ega_file_id,
            path_credentials_json=args.ega_credentials,
            executable=args.ega_executable,
            outdir=pcawg_file_id
        )
        bam = glob.glob(f"{pcawg_file_id}/{ega_file_id}/*.bam")[0]
    if result_download != 0:
        print(f"# download failed for pcawg_file_id={pcawg_file_id}, ega_file_id={ega_file_id}")
        return 1,0,0,0
        
    ### sort
    # sorted bam and all future steps in pipeline live up one dir from downloaded bam
    bam_sort = bam.replace(".bam", ".sort.bam")
    p = Path(bam_sort)
    bam_sort = "/".join([p.parts[0], p.parts[-1]])
    if args.cache and os.path.exists(bam_sort):
        result_sort = 0
    else:
        result_sort = samtools_sort(
            bam_in=bam,
            bam_out=bam_sort,
            threads=int(args.cpus_per_run)
        )
    if result_sort != 0:
        print(f"# sort failed for pcawg_file_id={pcawg_file_id}, bam={bam}")
        return 0,1,0,0
    ### bamtofastq
    r1_out = bam_sort.replace(".bam", ".r1.fq")
    r2_out = bam_sort.replace(".bam", ".r2.fq")
    if args.cache and os.path.exists(r1_out) and os.path.exists(r2_out):
        result_bamtofastq = 0
    else:
        result_bamtofastq = bamtofastq(
            bam_sort=bam_sort,
            r1_out=r1_out,
            r2_out=r2_out
        )
    if result_bamtofastq != 0:
        print(f"# bamtofastq failed for pcawg_file_id={pcawg_file_id}, bam_sort={bam_sort}")
        return 0,0,1,0
    ### star fusion
    breakpoint()
    outdir_star_fusion = os.path.join(os.path.dirname(bam_sort), args.outdir_star_fusion)
    if args.cache and os.path.isdir(outdir_star_fusion):
        result_star_fusion = 0
    else:
        result_star_fusion = star_fusion(
            fq1=r1_out,
            fq2=r2_out,
            outdir=outdir_star_fusion,
            star_fusion_conda=args.conda_env,
            cpus=int(args.cpus_per_run),
            genome_lib_dir=args.genome_lib_dir
        )
    if result_star_fusion != 0:
        print(f"# star fusion failed for pcawg_file_id={pcawg_file_id}, r1_out={r1_out},r2={r2_out}")
        return 0,0,0,1
    ### cleanup
    Path(bam).unlink(missing_ok=True)
    Path(bam_sort).unlink(missing_ok=True)
    Path(r1_out).unlink(missing_ok=True)
    Path(r2_out).unlink(missing_ok=True)
    return result_download, result_sort, result_bamtofastq, result_star_fusion

def main():
    validate_args()
    validate_executables()
    df = pd.read_csv(args.ega_table, sep="\t")
    # some (~54/577) of the RNA bams are missing from metadata file
    # we drop them because we cannot resolve the specimen type (e.g., tumor/normal) nor tissue 
    df.dropna(subset=['pcawg_file_id'],inplace=True)
    df.reset_index(drop=True, inplace=True)
    rows = [row for _, row in df.iterrows()]
    with Pool(processes=args.cpus_parallel) as pool:
        results = pool.map(process_row, rows)
    failed_rows = []
    for i,res_tuple in enumerate(results):
        # any will return True if any value is non-zero in the tuple
        if any(res_tuple):
            failed_rows.append(i)
    df_failed = df.iloc[failed_rows, :]
    df_failed.reset_index(drop=True,inplace=True)
    df_failed['pipeline_result_codes'] = [
        results[i] for i in failed_rows
    ]
    df_failed.to_csv(args.failed_rows,sep="\t",index=False)
    
    
        

if __name__ == "__main__":
    main()