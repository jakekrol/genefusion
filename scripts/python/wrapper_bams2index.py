#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys
import subprocess
import tempfile
import re
import multiprocessing as mp
import shutil
import time
import boto3

VALID_STEPS=range(7)

parser = argparse.ArgumentParser(description='Wrapper for building a GIGGLE index from discordant paired-end reads from BAMs')
parser.add_argument('-dd', '--discordant_dist', type=int, default=500)
parser.add_argument('-db', '--dirbam', type=str, required=True, help='Dir with BAM files')
parser.add_argument('-do', '--dirout', type=str, required=True, help='Dir to write output files')
parser.add_argument('-s', '--steps', type=str, default='all', help='Steps to run, comma separated')
parser.add_argument('-p', '--processes', type=int, default=1, help='Number of processes to use')
parser.add_argument('-dx', '--direxec', type=str, default='executables', help='Directory with executables' )
parser.add_argument('-bp', '--bam_pattern', type=str, default=None, help='BAM file pattern to match in dirbam')
parser.add_argument('--aws', action='store_true', help='Use AWS S3 for input BAM files')
parser.add_argument('--fileidmap', type=str, default='/data/jake/genefusion/data/meta/filename2fileid.tsv', help='File with mapping of filenames to file IDs')
parser.add_argument('-b', '--bedtools', type=str, default='/data/jake/bedtools.static.binary', help='Path to bedtools executable')
parser.add_argument('-ip', '--index_prefix', type=str, required=True, help='Prefix for index dir')
args = parser.parse_args()


def s3_download_file(bucket, src, dest):
    s3 = boto3.client('s3')
    s3.download_file(bucket, src, dest)

if not args.aws:
    assert os.path.exists(args.dirbam), f"Directory with BAM files does not exist: {args.dirbam}"
assert os.path.exists(args.dirout), f"Output directory does not exist: {args.dirout}"
assert os.path.exists(args.direxec), f"Directory with executables does not exist: {args.direxec}"
assert args.discordant_dist > 0, f"Discordant distance must be a positive integer, got: {args.discordant_dist}"

# parse steps
if not args.steps == 'all':
    # parse a range of steps e.g., 3- implies 3 to last step
    if "-" in args.steps:
        start = args.steps.split('-')[0]
        end = len(VALID_STEPS)
        args.steps = [int(i) for i in range(int(start), end)]
    # otherwise, parse a list of steps
    else: 
        args.steps = [int(i) for i in args.steps.split(',')]
    # validate
    if not all([i in VALID_STEPS for i in args.steps]):
        raise ValueError(f"Invalid steps: {args.steps}. Valid steps are: {VALID_STEPS}")
else:
    args.steps = list(VALID_STEPS)
print('Args:', args)
print(f"Running steps: {args.steps} with {args.processes} processes")
if args.aws:
    args.dirbam = 's3:/' + args.dirbam + '/'

# check executables
result = subprocess.run("gargs -h", shell=True,capture_output=True, text=True)
if result.returncode != 0:
    raise FileNotFoundError("gargs executable not found or not working. Please install gargs.")
if "4" in args.steps:
    if not os.path.exists(args.bedtools):
        raise FileNotFoundError(f"Bedtools executable not found: {args.bedtools}")

### inputs
if 0 in args.steps:
    # out subdirs
    subdirs = [
        "inputs",
        "exec",
        "logs",
        "bais",
        "xcord", # excord
        "xcord_cln", # clean
        "xcord_sort_zip" # sort and compress
        f"{args.index_prefix}_sort" # index
    ]
    print('Creating output subdirectories')
    for subdir in subdirs:
        os.makedirs(os.path.join(args.dirout, subdir), exist_ok=True)
    # executables
    print('Creating symlinks for executables')
    paths = [os.path.join(args.direxec, f) for f in os.listdir(args.direxec)]
    for src in paths:
        try:
            os.symlink(src, os.path.join(args.dirout, 'exec', os.path.basename(src)))
        except shutil.SameFileError:
            print(f"Symlink already exists for {src}, skipping.")
        except Exception as e:
            print(f"Error creating symlink for {src}: {e}")

    # BAMs
    print('Listing BAM files')
    if args.aws:
        print('Using AWS S3 for input BAM files:', args.dirbam)

        cmd = (
            f"aws s3 ls {args.dirbam} | awk '{{print $4}}' | sed '/^$/d' > {os.path.join(args.dirout, 'inputs', 'bams.list')}"
        )
        print('Running command to list BAM files:', cmd)
        subprocess.run(cmd, shell=True)
    else:
        sys.exit(f"Local BAM files are not supported yet. Please use --aws option to specify S3 bucket.")
    # filter BAMs by pattern if provided
    if args.bam_pattern:
        print(f"Filtering BAM files by pattern: {args.bam_pattern}")
        with open(os.path.join(args.dirout, 'inputs', 'bams.list'), 'r') as f:
            bams = [line.strip() for line in f if (args.bam_pattern in line.strip()) and line.strip().endswith('.bam')]
        with open(os.path.join(args.dirout, 'inputs', 'bams.list'), 'w') as f:
            for bam in bams:
                f.write(f"{bam}\n")
    # map file names to file IDs
    cmd = (
        f"./exec/pcawgbam2fileid.py "
        f"-f {args.fileidmap} "
        f"-i {os.path.join(args.dirout, 'inputs', 'bams.list')} "
        f"--suffix .excord.bed > {os.path.join(args.dirout, 'inputs', 'excord.input')}"
    )
    print('Running command:', cmd)
    subprocess.run(cmd, shell=True)
    # add dirout prefix to excord.input
    df = pd.read_csv(os.path.join(args.dirout, 'inputs', 'excord.input'), sep='\t', header=None)
    df[1] = df[1].apply(lambda x: os.path.join(args.dirout, 'xcord', x))
    df.to_csv(os.path.join(args.dirout, 'inputs', 'excord.input'), sep='\t', header=False, index=False)
    # make excord_cln.input
    df = pd.read_csv(os.path.join(args.dirout, 'inputs', 'excord.input'), sep='\t', header=None,usecols=[1])
    df.columns = ['excord']
    # for some reason excord outputs have .gz suffix, despite not being compressed
    # we'll keep the extension 
    df['excord'] = df['excord'].apply(lambda x: x + '.gz' if not x.endswith('.gz') else x)
    df['excord_cln'] = df['excord'].apply(
        lambda x: os.path.join(
            f'{args.dirout}/xcord_cln',
            os.path.basename(x)
        )
    )
    df.to_csv(os.path.join(args.dirout, 'inputs', 'excord_cln.input'), sep='\t', header=False, index=False)
    # make bedtools input
    df = pd.read_csv(os.path.join(args.dirout, 'inputs', 'excord_cln.input'), sep='\t', header=None, usecols=[1])
    df.columns = ['excord_cln']
    df['bedtools'] = df['excord_cln'].apply(
        lambda x: os.path.join(
            f'{args.dirout}/xcord_sort_zip',
            os.path.basename(x)
        )
    )
    df.to_csv(os.path.join(args.dirout, 'inputs', 'bedtools.input'), sep='\t', header=False, index=False)
if 1 in args.steps:
    # download BAIs
    print('Downloading BAI files')
    cmd = [
        "aws", "s3", "cp", args.dirbam, os.path.join(args.dirout, 'bais'),
        "--recursive", "--exclude", "*", "--include", "*.bai"
    ]
    print('Running command:', ' '.join(cmd))
    t = time.time()
    subprocess.run(cmd, check=True)
    print(f"Finished downloading BAI files in {time.time() - t:.2f} seconds")

if 2 in args.steps:
    # excord
    if args.aws:
        print('Running excord')
        cmd = (
            f"( cd {os.path.join(args.dirout, 'bais')} && "
            f"gargs -p {args.processes} "
            f"--log={os.path.join(args.dirout, 'logs', 'excord.log')} "
            f"-o 'aws s3 cp {args.dirbam}{{0}} - | excord --discordantdistance {str(args.discordant_dist)} /dev/stdin > {{1}}' "
            f"< {os.path.join(args.dirout, 'inputs', 'excord.input')} "
            f")"
        )
        print('Running command:', cmd)
        t= time.time()
        subprocess.run(cmd, shell=True)
        print(f"Finished excord in {time.time() - t:.2f} seconds")

if 3 in args.steps:
    # clean excord
    print('Cleaning excord files')
    cmd = (
        f"gargs -p {args.processes} "
        f"--log={os.path.join(args.dirout, 'logs', 'excord_clean.log')} "
        f"-o './exec/cln_excord.sh {{0}} {{1}}' "
        f"< {os.path.join(args.dirout, 'inputs', 'excord_cln.input')}"
    )
    print('Running command:', cmd)
    t = time.time()
    subprocess.run(cmd, shell=True)
    print(f"Finished cleaning excord files in {time.time() - t:.2f} seconds")   
if 4 in args.steps:
    # sort and compress excord
    print('Sorting and compressing excord files')
    cmd = (
        f"gargs -p {args.processes} "
        f"--log={os.path.join(args.dirout, 'logs', 'excord_sort_zip.log')} "
        f"-o '{args.bedtools} sort -i {{0}} | bgzip -c > {{1}}' "
        f"< {os.path.join(args.dirout, 'inputs', 'bedtools.input')}"
    )
    print('Running command:', cmd)
    t = time.time()
    subprocess.run(cmd, shell=True)
    print(f"Finished sorting and compressing excord files in {time.time() - t:.2f} seconds")
if 5 in args.steps:
    # copy sorted and compressed excord files to index dir
    print('Moving sorted and compressed excord files to index dir')
    files = os.listdir(os.path.join(args.dirout, 'xcord_sort_zip'))
    for f in files:
        src = os.path.join(args.dirout, 'xcord_sort_zip', f)
        dest = os.path.join(args.dirout, f"{args.index_prefix}_sort", f)
        shutil.copy(src, dest)



    