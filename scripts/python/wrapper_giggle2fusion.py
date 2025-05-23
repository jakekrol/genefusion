#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys
import subprocess
import tempfile
import re
import multiprocessing
import shutil
import time

# future refactoring:
# i) don't rely on filename to parse info. use a commented header instead
# ii) write each step to distinct log file

TEMPLATE_GIGGLE_SEARCH="/data/jake/genefusion/data/giggle_search.template"
FILEID2SAMPLETYPE="/data/jake/genefusion/data/meta/fileid2sampletype.tsv"
SAMPLECOLIDX = 15
SPECIMENCOLIDX = 16

parser = argparse.ArgumentParser(description="Run giggle2fusion pipeline")
parser.add_argument('-d', '--base_dir', type=str, required=True)
parser.add_argument('-g', '--giggle_index', type=str, required=True)
parser.add_argument('-s', '--step', type=int, default=0, help='Step to start from')
parser.add_argument('-x', '--dir_executables', type=str, default='/data/jake/genefusion/executables', help='Directory containing executables')  
parser.add_argument('-p', '--processes', type=int, default=60, help='Number of processes to use')
parser.add_argument('-f', '--freeze', action='store_true', help='Freeze: only run the specified step and exit')
args = parser.parse_args()
print("Running giggle2fusion pipeline from step: ", args.step)
time.sleep(3)


# def cln_sample_names(x):
#     # remove index prefix
#     x = os.path.basename(x)
#     x = re.sub("r\.excord\.bed\.gz", "", x)
#     return x

steps = {
    0: "inputs",
    1: "search",
    2: "clean",
    3: "swap",
    4: "intersect",
    5: "unswap",
    6: "unswap_cln",
    7: "unswap_specimen",
    8: "unswap_specimen_split",
    9: "migrate_specimen"
}

# scripts = [
#     "genefusion_giggle.sh",
#     "cln_excord.sh",
#     "swap_intervals.sh",
#     "intersect_swapped.sh",
#     "unswap_intervals.sh",
#     "cln_sample_name.py",
#     "add_specimentype.py",
#     "specimensplit_intersect.sh",
#     "migrate_specimen.py",
#     "count_fusions.sh",
#     "build_fusionpe_tbl2.py",
#     "count_samples_w1.py",
#     "burden_total.py"
# ]

### inputs
if args.step == 0:
    t = time.time()
    # make subdirs
    subdirs = [
        "inputs",
        "giggleout",
        "gigglecln",
        "giggleswap",
        "giggleinter",
        "giggleinter_unswap",
        "giggleinter_unswap_cln",
        "giggleinter_unswap_specimen",
        "giggleinter_unswap_specimen_split"
        "giggleinter_final_tumor",
        "giggleinter_final_normal"
    ]
    print('Making subdirs')
    for subdir in subdirs:
        os.makedirs(os.path.join(args.base_dir, subdir), exist_ok=True)
    # copy executables
    print(f'Copying executables from {args.dir_executables}')
    paths = [os.path.join(args.dir_executables, i) for i in os.listdir(args.dir_executables)]
    for src in paths:
        # copy to base_dir
        try:
            shutil.copy(src, args.base_dir)
            print(f"Copied {src} to {args.base_dir}")
        except shutil.SameFileError:
            print(f"File {src} already exists in {args.base_dir}")
        except Exception as e:
            print(f"An error occurred: {e}")
# make search input file
    print(f"Making search input file: {TEMPLATE_GIGGLE_SEARCH}")
    df_search = pd.read_csv(TEMPLATE_GIGGLE_SEARCH, sep="\t", header=None)
    df_search.columns = ["chrom", "start", "end", "name", "strand", "filename"]
    if '/' in args.giggle_index: # assume not in base_dir
        df_search.insert(0, "giggle_index", args.giggle_index)
        print (f"Using giggle index: {args.giggle_index}")
    else:
        df_search.insert(0, "giggle_index", os.path.join(args.base_dir,args.giggle_index))
        print (f"Using giggle index: {os.path.join(args.base_dir,args.giggle_index)}")
    fnames = df_search['filename'].tolist()
    df_search["filename"] = df_search["filename"].apply(lambda x: os.path.join(args.base_dir, "giggleout", x))
    df_search.to_csv(os.path.join(args.base_dir, "inputs", "search.input"), sep="\t", index=False, header=False)

    # make clean input file
    input_file = os.path.join(args.base_dir, "inputs", "clean.input")
    print(f"Making clean input file: {input_file}")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleout', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'gigglecln', fname)}\n")
    # make swap input file
    input_file = os.path.join(args.base_dir, "inputs", "swap.input")
    print(f"Making swap input file: {input_file}")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'gigglecln', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'giggleswap', fname)}\n")
    # make intersect input file
    input_file = os.path.join(args.base_dir, "inputs", "intersect.input")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleswap', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'giggleinter', fname)}\n")
    # make unswap input file
    input_file = os.path.join(args.base_dir, "inputs", "unswap.input")
    print(f"Making unswap input file: {input_file}")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleinter', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_unswap', fname)}\n")
    # make clean unswap input file
    input_file = os.path.join(args.base_dir, "inputs", "unswap_cln.input")
    print(f"Making unswap specimen input file: {input_file}")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_unswap', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_unswap_cln', fname)}\n")
    # make specimen input file
    input_file = os.path.join(args.base_dir, "inputs", "unswap_specimen.input")
    print(f"Making unswap specimen input file: {input_file}")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_unswap_cln', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_unswap_specimen', fname)}\n")
    # make specimen split input file
    input_file = os.path.join(args.base_dir, "inputs", "unswap_specimen_split.input")
    print(f"Making unswap specimen split input file: {input_file}")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_unswap_specimen', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_unswap_specimen_split')}\n") # constant output dir
    print(f"Finished making input files in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1

### run
if args.step == 1:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "search.input")), "Search input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleout")), "Giggle output directory not empty"
    cmd = [
    "gargs",
    "-p", f"{args.processes}",
    "--log=g.log",
    "-o", "./genefusion_giggle.sh -i {0} -g {4} -c {1} -l {2} -r {3} -s {5} -o {6}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "search.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    print(f"Finished running search in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1

if args.step == 2:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "clean.input")), "Clean input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "gigglecln")), "Giggle clean directory not empty"
    cmd = [
    "gargs",
    "-p", f"{args.processes}",
    "--log=g.log",
    "-o", "./cln_excord.sh {0} {1}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "clean.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    print(f"Finished running clean in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1
if args.step == 3:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "swap.input")), "Swap input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleswap")), "Giggle swap directory not empty"
    cmd = [
    "gargs",
    "-p", f"{args.processes}",
    "--log=g.log",
    "-o", "./swap_intervals.sh {0} {1}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "swap.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t= time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    print(f"Finished running swap in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1
if args.step == 4:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "intersect.input")), "Intersect input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleinter")), "Giggle intersect directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", "./intersect_swapped.sh {0} {1}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "intersect.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1
if args.step == 5:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "unswap.input")), "Unswap input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleinter_unswap")), "Giggle unswap directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", "./unswap_intervals.sh {0} {1}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "unswap.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    print(f"Finished running unswap in {time.time() - t:.2f} seconds")
if args.step == 6:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "unswap_cln.input")), "Unswap clean input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleinter_unswap_cln")), "Giggle unswap clean directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", f"./cln_sample_name.py -i {{0}} -o {{1}} -s {SAMPLECOLIDX}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "unswap_cln.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    print(f"Finished running unswap clean in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1
if args.step == 7:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "unswap_specimen.input")), "Unswap specimen input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleinter_unswap_specimen")), "Giggle unswap specimen directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", f"./add_specimentype.py -i {{0}} -o {{1}} -s {SAMPLECOLIDX} -l {FILEID2SAMPLETYPE}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "unswap_specimen.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    print(f"Finished running unswap specimen in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1
if args.step == 8:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "unswap_specimen_split.input")), "Unswap specimen split input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleinter_unswap_specimen_split")), "Giggle unswap specimen split directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", f"./specimensplit_intersect.sh -i {{0}} -o {{1}} -s {SPECIMENCOLIDX}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "unswap_specimen_split.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
        subprocess.run(cmd)
    print(f"Finished running unswap specimen split in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1
if args.step == 9:
    assert os.path.isdir(os.path.join(args.base_dir, "giggleinter_final_tumor")), "giggleinter_final_tumor directory not found"
    assert os.path.isdir(os.path.join(args.base_dir, "giggleinter_final_normal")), "giggleinter_final_normal directory not found"
    cmd = [
        "./migrate_specimen.py",
        "-i", os.path.join(args.base_dir, "giggleinter_unswap_specimen_split"),
        "-t", os.path.join(args.base_dir, "giggleinter_final_tumor"),
        "-n", os.path.join(args.base_dir, "giggleinter_final_normal"),
        "-p", f"{args.processes}"
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running migrate specimen in {time.time() - t:.2f} seconds")
    if args.freeze:
        print(f"Freezing at step {args.step}")
        sys.exit(0)
    args.step +=1

