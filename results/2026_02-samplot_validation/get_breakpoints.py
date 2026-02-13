#!/usr/bin/env python3
import argparse
import pandas as pd
import os,sys
import subprocess
import tempfile

parser=argparse.ArgumentParser(description='Get breakpoints of fusions using giggle')
parser.add_argument('-i','--input',help='Input fusion file',required=True)
parser.add_argument('--index_root', default='/data/jake/stix-pcawg-dna')
parser.add_argument('-o', '--output', default='fusion_breakpoints.tsv')
parser.add_argument('-t', '--tmpdir', default='/tmp')
parser.add_argument('--swap_script', 
                    default='../../scripts/shell/swap_intervals.sh')
parser.add_argument('--intersect_swap_script',
                    default='../../scripts/shell/intersect_swapped.sh')
parser.add_argument('-b', '--gene_bedfile',
                    default='../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed'
)
args=parser.parse_args()
# globally set temp dir
tempfile.tempdir = args.tmpdir

df_in = pd.read_csv(args.input, sep='\t')
df_in['left_breakpoint'] = pd.NA
df_in['right_breakpoint'] = pd.NA
for i,row in df_in.iterrows():

    # get info for giggle search
    shard = row['shard']
    generight = row['right']
    dirindex= os.path.join(args.index_root, shard)
    left_chrom=row['left_gene_chrom']
    left_start=row['left_gene_start']
    left_end=row['left_gene_end']
    region_str=f'{left_chrom}:{left_start}-{left_end}'

    # run giggle search
    tmp_giggle = tempfile.NamedTemporaryFile(mode='w', delete=False)
    cmd_giggle = f'( cd {args.index_root} && giggle search -i {dirindex} ' \
        f'-r {region_str} -v > {tmp_giggle.name} )'
    subprocess.run(cmd_giggle, shell=True)
    # run interval column swap
    tmp_swap_out = tempfile.NamedTemporaryFile(mode='w', delete=False)
    cmd_swap= f'{args.swap_script} {tmp_giggle.name} {tmp_swap_out.name}'
    subprocess.run(cmd_swap, shell=True)
    # run intersect swapped
    tmp_intersect_out = tempfile.NamedTemporaryFile(mode='w', delete=False)
    cmd_intersect = f'{args.intersect_swap_script} {tmp_swap_out.name} ' \
        f'{args.gene_bedfile} {tmp_intersect_out.name}'
    subprocess.run(cmd_intersect, shell=True)
    # run right gene filter
    tmp_filtered_out = tempfile.NamedTemporaryFile(mode='w', delete=False)
    cmd_filter = f'grep -w {generight} {tmp_intersect_out.name} > {tmp_filtered_out.name}'
    subprocess.run(cmd_filter, shell=True)
    # get breakpoints
    cmd_min_left = f'cut -f 11 {tmp_filtered_out.name} | sort -n | head -n 1'
    cmd_max_right = f'cut -f 8 {tmp_filtered_out.name} | sort -nr | head -n 1'
    left_bp = subprocess.check_output(cmd_min_left, shell=True).decode('utf-8').strip()
    right_bp = subprocess.check_output(cmd_max_right, shell=True).decode('utf-8').strip()
    df_in.at[i,'left_breakpoint'] = left_bp
    df_in.at[i,'right_breakpoint'] = right_bp
    # cleanup temp files
    os.remove(tmp_giggle.name)
    os.remove(tmp_swap_out.name)
    os.remove(tmp_intersect_out.name)
    os.remove(tmp_filtered_out.name)

# write
df_in.to_csv(args.output, sep='\t', index=False)
