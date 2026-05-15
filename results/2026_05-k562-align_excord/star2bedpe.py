#!/usr/bin/env python3
from polymerization.star import *
from polymerization.utils import *
from polymerization.giggle2fusion import *
import argparse
import os
import tempfile
import time
parser=argparse.ArgumentParser(description='Convert STAR Chimeric.out.junction format to sorted BEDPE format')
parser.add_argument('--input', default='k562-rna-fusion-Chimeric.out.junction', help='Path to STAR Chimeric.out.junction file')
parser.add_argument('--output', default='k562-rna-fusion.bedpe.gz', help='Output BEDPE file path (should end with .gz if --bgzip is set)')
parser.add_argument('--has_header', default=1, help='Input file has a header line')
parser.add_argument('--bgzip', default=1, help='Compress output with bgzip')
parser.add_argument('--tmp_dir', default='/data/jake/tmp', help='Temporary directory for sorting')
args=parser.parse_args()
for arg in vars(args):
    print(f"{arg}: {getattr(args, arg)}")

with tempfile.NamedTemporaryFile() as tmp_bedpe:
    tmp_bedpe_path = tmp_bedpe.name
    t_0=time.time()
    chimeric2bedpe(args.input, tmp_bedpe_path, has_header=args.has_header, bgzip=False)
    with tempfile.NamedTemporaryFile() as tmp_sorted_bedpe:
        tmp_sorted_bedpe_path = tmp_sorted_bedpe.name + '.gz' if args.bgzip else tmp_sorted_bedpe.name
        sort_bed(tmp_bedpe_path, tmp_sorted_bedpe_path, bgzip=args.bgzip)
        print(f'# time to convert STAR chimeric to BEDPE and sort: {time.time() - t_0:.2f} seconds')
        clean_excord(tmp_sorted_bedpe_path, args.output, bgzip=args.bgzip)


