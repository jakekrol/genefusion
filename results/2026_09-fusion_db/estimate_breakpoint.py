#!/usr/bin/env python3

import argparse
import glob
import os
from multiprocessing import Pool
import pandas as pd
from polymerization.analysis import intersect_file2breakpoints
from polymerization.datasets import get_pcawg_data_types
from polymerization.io import read_bed, read_g2f_intersect

parser = argparse.ArgumentParser()
parser.add_argument("--dir_g2f", default="../2026_06-g2f-all_gene_pairs/g2f_out")
parser.add_argument("--indir", default='top_fusions')
parser.add_argument("--infile_suffix", default='.annotate.repeat_filt.tsv')
parser.add_argument("--outdir",default="top_fusions")
parser.add_argument("--outfile_suffix", default=".annotate.repeat_filt.breakpoint.tsv")
parser.add_argument("--intersect_file_suffix", default=".giggle.clean.swap.intersect.bed.gz")
parser.add_argument("--bed", default = "../2026_09-living_bed/grch37.genes.sort.bed")
parser.add_argument("--cpus", default = 30)
args = parser.parse_args()
args.cpus = int(args.cpus)
MODALITY='rna'

def get_pcawg_tumor_tissues():
    data = get_pcawg_data_types()
    tissues = set()
    for tissue, _ in data.items():
        # if 'tumor' in data[tissue]['modality']['dna']['specimen']:
        #     tissues.add(tissue)
        if 'tumor' in data[tissue]['modality']['rna']['specimen']:
            tissues.add(tissue)
    return tissues

def validate_args():
    assert os.path.exists(args.dir_g2f)
    assert os.path.exists(args.indir)
    assert args.cpus > 0

def find_intersect_file(gene, tissue, directory=args.dir_g2f, suffix=args.intersect_file_suffix):
    directory = os.path.join(
        args.dir_g2f,
        f"{tissue}_tumor_{MODALITY}"
    )
    intersect_file = os.path.join(
        directory,
        f"{gene}{suffix}"
    )
    if os.path.exists(intersect_file):
        return intersect_file
    else:
        return None

def estimate_breakpoint(in_filepath: str, out_filepath: str, tumor_tissues: list):
    tissue = os.path.basename(in_filepath).split('-')[3]
    if tissue not in tumor_tissues:
        return None
    df = pd.read_csv(in_filepath,sep='\t')
    df['breakpoint_left'] = pd.NA
    df['breakpoint_right'] = pd.NA
    for i,row in df.iterrows():
        gene_left = row['gene_left']
        gene_right = row['gene_right']
        intersect_file = find_intersect_file(gene_left,tissue)
        if not intersect_file:
            continue
        df_bp = intersect_file2breakpoints(
            path=intersect_file,
            gene_right_query=gene_right,
            group_by_sample=True,
            bgzip=True
        )
        # intersect_file2breakpoints failed due to timeout or empty file most likely
        if type(df_bp) != type(pd.DataFrame()):
            print(f"# failed to estimate breakpoint for {gene_left}--{gene_right}")
            continue
        left_chrom = df_bp.loc[0,'chromosome_left']
        right_chrom = df_bp.loc[0,'chromosome_right']
        bp_left = int(df_bp['start_left'].mean())
        bp_right = int(df_bp['end_right'].mean())
        df.at[i, 'breakpoint_left'] = f"{left_chrom}:{bp_left}"
        df.at[i, 'breakpoint_right'] = f'{right_chrom}:{bp_right}'
    df.to_csv(out_filepath,sep='\t',index=False)

def main():
    tumor_tissues = {"blood", "kidney", "liver", "ovary"}
    os.makedirs(args.outdir,exist_ok=True)
    infiles = glob.glob(f"{args.indir}/*{args.infile_suffix}")
    infiles.sort()
    queries = []
    for f in infiles:
        for tissue in tumor_tissues:
            if tissue in f:
                outfile = f.replace(args.infile_suffix,args.outfile_suffix)
                queries.append((f,outfile,tumor_tissues))
    cpus = min(args.cpus, len(queries))
    print(f"# cpus: {cpus}")

    with Pool(cpus) as pool:
        try:
            results = pool.starmap(estimate_breakpoint, queries)
        except KeyboardInterrupt:
            pool.terminate()
            raise

if __name__ == "__main__":
    main()