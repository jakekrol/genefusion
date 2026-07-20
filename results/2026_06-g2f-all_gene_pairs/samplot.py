#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import multiprocessing as mp
import subprocess
import shlex

parser = argparse.ArgumentParser(description="")
parser.add_argument("--bams", default='bams', help="Input directory path")
parser.add_argument("--outdir", default='samplots', type=str, help="Output directory path")
parser.add_argument("--tissue", type=str, help="Tissue name")
parser.add_argument("--table", type=str, help="Input table with fusion information")
parser.add_argument("--gene_bed", default='/data/jake/jkbiolib/jkbiolib/data/gencode.v19.annotation.gtf.gene.bed.sorted.gz', type=str, help="Gene annotation bed file")
parser.add_argument('--threads', type=int, default=8, help='Number of threads to use for multiprocessing')
args = parser.parse_args()
os.makedirs(args.outdir, exist_ok=True)

def make_annotation_bed(annotation_bed_path, annotation_bed_out, gene_left, gene_right):
    # Simple, short pipeline: zcat -> grep -wE -> bgzip -> tabix
    pattern = f"{gene_left}|{gene_right}"
    cmd = (
        f"zcat {shlex.quote(annotation_bed_path)} | grep -wE {shlex.quote(pattern)} | bgzip -c > {shlex.quote(annotation_bed_out)} && "
        f"tabix -p bed {shlex.quote(annotation_bed_out)}"
    )
    print("# running command: {}".format(cmd))
    subprocess.run(cmd, shell=True, check=True)
    return annotation_bed_out


def samplot(bam, outfile, name, annotation_bed, chromosome_list, start_list, end_list):
    if (len(chromosome_list) == 2) and (len(start_list) == 2) and (len(end_list) == 2):
        chr1 = chromosome_list[0]
        start1 = start_list[0]
        end1 = end_list[0]
        chr2 = chromosome_list[1]
        start2 = start_list[1]
        end2 = end_list[1]
        cmd = [
            'samplot', 'plot',
            '-b', bam,
            '-o', outfile,
            '-n', name,
            '-A', annotation_bed,
            '-c', str(chr1),
            '-s', str(start1),
            '-e', str(end1),
            '-c', str(chr2),
            '-s', str(start2),
            '-e', str(end2),
            '-t', 'BND'
        ]
        print("# running command: {}".format(' '.join(cmd)))
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = result.returncode
        return retcode
    else:
        chrom = chromosome_list[0]
        start = start_list[0]
        end = end_list[0]
        cmd = [
            'samplot', 'plot',
            '-b', bam,
            '-o', outfile,
            '-n', name,
            '-A', annotation_bed,
            '-c', str(chrom),
            '-s', str(start),
            '-e', str(end)
        ]
        print("# running command: {}".format(' '.join(cmd)))
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = result.returncode
        return retcode
def prep_bam(bam_in, bam_out):
    # sort
    cmd_sort = ['samtools', 'sort', '-o', bam_out, bam_in]
    print("# running command: {}".format(' '.join(cmd_sort)))
    result_sort = subprocess.run(cmd_sort, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # index
    cmd_index = ['samtools', 'index', '-o', f"{bam_out}.bai", bam_out]
    print("# running command: {}".format(' '.join(cmd_index)))
    result_index = subprocess.run(cmd_index, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    
    
df = pd.read_csv(args.table, sep='\t')
for i, row in df.iterrows():
    fusion = row['fusion']
    gene_left = row['gene_left']
    gene_right = row['gene_right']
    sample = row['sample']
    specimen = row['specimen']
    bam = row['region_bam']
    bam = os.path.basename(bam)
    bam = os.path.join(args.bams, args.tissue, bam)
    outfile = os.path.join(args.outdir, f"{fusion}.{sample}.png")
    name = f"{fusion};{sample};{specimen}"
    left_chromosome = row['left_chromosome']
    right_chromosome = row['right_chromosome']
    left_bp = row['left_breakpoint']
    right_bp = row['right_breakpoint']
    if left_chromosome == right_chromosome:
        chromosome_list = [left_chromosome]
        start_list = [left_bp]
        end_list = [right_bp]
    else:
        chromosome_list = [left_chromosome, right_chromosome]
        start_list = [left_bp, left_bp]
        end_list = [right_bp, right_bp]
    annotation_bed_subset = os.path.join(args.outdir, f"{fusion}.bed.gz")
    bam_sort = bam.replace('.bam', '.sorted.bam')
    try:
        annotation_bed_out = make_annotation_bed(args.gene_bed, annotation_bed_subset, gene_left, gene_right)
        result = prep_bam(bam, bam_sort)
        result = samplot(bam_sort, outfile, name, annotation_bed_out, chromosome_list, start_list, end_list)
    except Exception as e:
        print("# Error processing fusion {}: {}".format(fusion, str(e)))

