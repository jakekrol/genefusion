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
# ii) add a logfile

TEMPLATE_GIGGLE_SEARCH="/data/jake/genefusion/data/giggle_search.template"
FILEID2SAMPLETYPE="/data/jake/genefusion/data/meta/fileid2sampletype.tsv"
SAMPLECOLIDX = 15
SPECIMENCOLIDX = 16
VALID_STEPS = range(18)


parser = argparse.ArgumentParser(description="Run giggle2fusion pipeline")
parser.add_argument('-d', '--base_dir', type=str, required=True)
parser.add_argument('-g', '--giggle_index', type=str, required=True)
parser.add_argument('-s', '--steps', type=str, default='all', help='Steps to run, comma separated')
parser.add_argument('-x', '--dir_executables', type=str, default='/data/jake/genefusion/executables', help='Directory containing executables')  
parser.add_argument('-p', '--processes', type=int, default=60, help='Number of processes to use')
args = parser.parse_args()
print('Running wrapper_giggle2fusion.py with arguments:')
print(args)
time.sleep(3)

# parse steps
if not args.steps == 'all':
    args.steps = [int(i) for i in args.steps.split(',')]
    if not all([i in VALID_STEPS for i in args.steps]):
        raise ValueError(f"Invalid steps: {args.steps}. Valid steps are: {VALID_STEPS}")
else:
    args.steps = list(VALID_STEPS)
print(f"Steps to run: {args.steps}")
time.sleep(3)


# def cln_sample_names(x):
#     # remove index prefix
#     x = os.path.basename(x)
#     x = re.sub("r\.excord\.bed\.gz", "", x)
#     return x


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
if 0 in args.steps:
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
        "giggleinter_unswap_specimen_split",
        "giggleinter_final_tumor",
        "giggleinter_final_normal",
        "pop_tumor_fusion_counts",
        "pop_tumor_fusion_sample_counts"
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
            os.symlink(src, args.base_dir)
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
    # make count fusions input file
    input_file = os.path.join(args.base_dir, "inputs", "count_fusions_tumor.input")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_final_tumor', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'pop_tumor_fusion_counts', fname)}\n")
    # make distinct sample counts input file
    input_file = os.path.join(args.base_dir, "inputs", "distinct_sample_counts.input")
    with open(input_file, 'w') as f:
        for fname in fnames:
            f.write(f"{os.path.join(args.base_dir, 'giggleinter_final_tumor', fname)}\t")
            f.write(f"{os.path.join(args.base_dir, 'pop_tumor_fusion_sample_counts', fname)}\n")

### run
if 1 in args.steps:
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
    print(f"Finished running search in {time.time() - t:.2f} seconds")

if 2 in args.steps:
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
    print(f"Finished running clean in {time.time() - t:.2f} seconds")

if 3 in args.steps:
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
    print(f"Finished running swap in {time.time() - t:.2f} seconds")

if 4 in args.steps:
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

if 5 in args.steps:
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
    print(f"Finished running unswap in {time.time() - t:.2f} seconds")
if 6 in args.steps:
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
    print(f"Finished running unswap clean in {time.time() - t:.2f} seconds")
if 7 in args.steps:
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
    print(f"Finished running unswap specimen in {time.time() - t:.2f} seconds")

if 8 in args.steps:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "unswap_specimen_split.input")), "Unswap specimen split input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "giggleinter_unswap_specimen_split")), "Giggle unswap specimen split directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", f"./specimensplit_intersect.sh {{0}} {{1}}"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "unswap_specimen_split.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
    print(f"Finished running unswap specimen split in {time.time() - t:.2f} seconds")

if 9 in args.steps:
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

if 10 in args.steps:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "count_fusions_tumor.input")), "Count fusions input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "pop_tumor_fusion_counts")), "Pop tumor fusion counts directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", f"./count_fusions.py -i {{0}} -o {{1}} -z -r 4"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "count_fusions_tumor.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
    print(f"Finished running count fusions in {time.time() - t:.2f} seconds")

if 11 in args.steps:
    assert os.path.isdir(os.path.join(args.base_dir, "pop_tumor_fusion_counts")), "pop_tumor_fusion_counts directory not found"
    assert not os.path.exists(os.path.join(args.base_dir, "pop_tumor_fusions.tsv")), "pop_tumor_fusions.tsv file already exists"
    cmd = [
        "./agg_pe_counts.py",
        "-i", os.path.join(args.base_dir, "pop_tumor_fusion_counts"),
        "-o", os.path.join(args.base_dir, "pop_tumor_fusions.tsv")
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running agg_pe_counts in {time.time() - t:.2f} seconds")

if 12 in args.steps:
    assert os.path.exists(os.path.join(args.base_dir, "inputs", "distinct_sample_counts.input")), "Distinct sample counts input file not found"
    assert not os.listdir(os.path.join(args.base_dir, "pop_tumor_fusion_sample_counts")), "Pop tumor fusion sample counts directory not empty"
    cmd = [
        "gargs",
        "-p", f"{args.processes}",
        "--log=g.log",
        "-o", f"./distinct_sample_counts.py -i {{0}} -o {{1}} -r 4 -s 15 -z"
    ]
    input_file = os.path.join(args.base_dir, "inputs", "distinct_sample_counts.input")
    print(f"Running '{' '.join(cmd)}' with input file: {input_file}")
    t = time.time()
    with open(input_file, 'r') as infile:
        subprocess.run(cmd, stdin=infile)
    print(f"Finished running distinct sample counts in {time.time() - t:.2f} seconds")
    # aggregate
    cmd = [
        "for",
        "file",
        "in",
        os.path.join(args.base_dir, "pop_tumor_fusion_sample_counts", "*"),
        ";",
        "do",
        "tail",
        "-n", "+2",
        "--quiet",
        "$file",
        ">>",
        os.path.join(args.base_dir, "pop_tumor_fusion_sample_counts.tsv"),
        ";",
        "done"
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(" ".join(cmd), shell=True)
    print(f"Finished running distinct sample counts aggregation in {time.time() - t:.2f} seconds")
    # add header to the aggregated file
    cmd = [
        "sed",
        "-i",
        "1ileft\tright\tsample_count",
        os.path.join(args.base_dir, "pop_tumor_fusion_sample_counts.tsv")
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished adding header in {time.time() - t:.2f} seconds")
if 13 in args.steps:
    # ./burden_total.py -i pop_tumor_fusions.tsv -o burden_total_tumor.tsv
    assert os.path.exists(os.path.join(args.base_dir, "pop_tumor_fusions.tsv")), "pop_tumor_fusions.tsv file not found"
    assert not os.path.exists(os.path.join(args.base_dir, "burden_total_tumor.tsv")), "burden_total_tumor.tsv file already exists"
    cmd = [
        "./burden_total.py",
        "-i", os.path.join(args.base_dir, "pop_tumor_fusions.tsv"),
        "-o", os.path.join(args.base_dir, "burden_total_tumor.tsv"),
        "--header"
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running burden_total in {time.time() - t:.2f} seconds")
if 14 in args.steps:
    assert os.path.exists(os.path.join(args.base_dir, "pop_tumor_fusions.tsv")), "pop_tumor_fusions.tsv file not found"
    assert os.path.exists(os.path.join(args.base_dir, "pop_tumor_fusion_sample_counts.tsv")), "pop_tumor_fusion_sample_counts.tsv file not found"
    cmd = [
        "./join.py",
        "-x", os.path.join(args.base_dir, "pop_tumor_fusions.tsv"),
        "-y", os.path.join(args.base_dir, "pop_tumor_fusion_sample_counts.tsv"),
        "-t", "left",
        "-o", os.path.join(args.base_dir, "pop_tumor_fusions_pe_and_sample.tsv"),
        "-k", "left,right"
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running join in {time.time() - t:.2f} seconds")
if 15 in args.steps:
    assert os.path.exists(os.path.join(args.base_dir, "pop_tumor_fusions_pe_and_sample.tsv")), "pop_tumor_fusions_pe_and_sample.tsv file not found"
    assert os.path.exists(os.path.join(args.base_dir, "burden_total_tumor.tsv")), "burden_total_tumor.tsv file not found"
    # left gene burden
    cmd = [
        "./add_burden_col.py",
        "-f", os.path.join(args.base_dir, "pop_tumor_fusions_pe_and_sample.tsv"),
        "-b", os.path.join(args.base_dir, "burden_total_tumor.tsv"),
        "-o", os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden_tmp.tsv"),
        "-k1", "0",
        "-k2", "0",
        "-n", "burden_total_left",
        "-hf",
        "-hb"
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running add_burden_col left in {time.time() - t:.2f} seconds")
    # right gene burden
    cmd = [
        "./add_burden_col.py",
        "-f", os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden_tmp.tsv"),
        "-b", os.path.join(args.base_dir, "burden_total_tumor.tsv"),
        "-o", os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden.tsv"),
        "-k1", "1",
        "-k2", "0",
        "-n", "burden_total_right",
        "-hf",
        "-hb"
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running add_burden_col right in {time.time() - t:.2f} seconds")
    # remove tmp file
    os.remove(os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden_tmp.tsv"))
    print(f"Removed temporary file: {os.path.join(args.base_dir, 'pop_tumor_fusions_pe_sample_burden_tmp.tsv')}")
if 16 in args.steps:
    assert os.path.exists(os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden.tsv")), "pop_tumor_fusions_pe_sample_burden.tsv file not found"
    cmd = [
        "./add_sample_density.py",
        "-i", os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden.tsv"),
        "-o", os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden_density.tsv")
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running add_sample_density in {time.time() - t:.2f} seconds")
if 17 in args.steps:
    assert os.path.exists(os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden_density.tsv")), "pop_tumor_fusions_pe_sample_burden_density.tsv file not found"
    cmd = [
        "./add_burden_product.py",
        "-i", os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden_density.tsv"),
        "-o", os.path.join(args.base_dir, "pop_tumor_fusions_pe_sample_burden_density_product.tsv")
    ]
    print(f"Running '{' '.join(cmd)}'")
    t = time.time()
    subprocess.run(cmd)
    print(f"Finished running add_burden_product in {time.time() - t:.2f} seconds")

