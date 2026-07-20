#!/usr/bin/env python3

import subprocess
import os
import argparse
import glob
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description='Merge BAM files')
parser.add_argument('-i', '--input_dir', default='bams/blood', help='Input BAMs')
parser.add_argument('-o', '--output_dir', default='bams_merge/blood', help='Output BAMs')
parser.add_argument('--threads', type=int, default=8, help='Number of threads to use')
args = parser.parse_args()
os.makedirs(args.output_dir, exist_ok=True)

bams = glob.glob(f"{args.input_dir}/*.bam")
bams = [os.path.basename(bam) for bam in bams]
fusions = set()
for bam in bams:
    fusion = bam.split('.')[:-3][0]
    fusions.add(fusion)
# Merge BAM files for each fusion
for fusion in fusions:
    print("# merging BAM files for fusion:", fusion)
    files_to_merge = [bam for bam in bams if bam.startswith(fusion)]
    n_files = len(files_to_merge)
    with tempfile.TemporaryDirectory(prefix="sort_bams_") as tmp:
        tmpdir = Path(tmp)
        # Sort BAM files before merging
        print(f"# sorting {n_files} BAM filesfor merge")
        sorted_bams = []
        for bam in files_to_merge:
            input_bam = Path(args.input_dir) / bam
            sorted_bam = tmpdir / f"{bam}.sorted.bam"
            sort_cmd = [
                'samtools', 'sort',
                '-@', str(args.threads),
                '-o', str(sorted_bam),
                str(input_bam)
            ]
            subprocess.run(sort_cmd, check=True)
            sorted_bams.append(str(sorted_bam))
        # Merge sorted BAM files
        print(f"# merging sorted BAM files to {args.output_dir}/{fusion}.bam")
        cmd = [
            'samtools', 'merge',
            '-@', str(args.threads),
            '--write-index',
            '-O', 'BAM',
            '-o', f"{args.output_dir}/{fusion}.bam",
        ] + sorted_bams
        subprocess.run(cmd, check=True)