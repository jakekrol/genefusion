#!/usr/bin/env python3

from genefusion.genefusion import get_sample_wise_fusions
import os, sys

# example 
# gf-samplefusions /data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/fusions_25_01_05/10.neg.A1CF.52559169.52645435.fusion /data/jake/genefusion/scratch/2025-01-17-sample_wise_fusions

# split fusions sample wise to file
args = sys.argv

genefusionfile = args[1]
outdir = args[2]
# get gene from filename
gene = os.path.basename(genefusionfile).split('.')[2:-3][0]
# will append to outdir/sample.fusion

def main():
    print("genefusionfile:", genefusionfile, "outdir: ", outdir, "gene: ", gene)
    output_file = os.path.join(outdir, f"{gene}.fusion")
    
    get_sample_wise_fusions(
        genefusionfile,
        gene,
        outdir,
        append=True
    )

if __name__ == "__main__":
    main()
