#!/usr/bin/env python

# LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH ./add_coords2stix_summary.py 
import os,sys
import pandas as pd
from genefusion.genefusion import left_gene
import argparse
import time

t_0=time.time()
parser=argparse.ArgumentParser(description='Add gene coordinates to STIX summary file.')
parser.add_argument('-i', '--input', help='Input STIX summary file', required=True)
parser.add_argument('-o', '--output', help='Output STIX summary file with gene coordinates', required=True)
parser.add_argument('-b', '--bedfile', help='Gene BED file with promoter padding', default='../2025_04-gene_bedfile_cln/grch37.genes.promoter_pad.bed')
args = parser.parse_args()

# validation
assert os.path.isfile(args.input), f"Input file {args.input} does not exist."
assert os.path.isfile(args.bedfile), f"BED file {args.bedfile} does not exist."
assert args.input != args.output, "Input and output files must be different."

df_bed=pd.read_csv(args.bedfile,sep='\t',header=None,names=['chrom','start','end','gene','strand'])
# parse gene names from filename
df_f=pd.read_csv(args.input,sep='\t')
gene1=args.input.split('.stix.')[1].split('--')[0]
gene2=args.input.split('.stix.')[1].split('--')[1].split('.out')[0]
# setup dataframe for left_gene function
df_genes=pd.DataFrame({'x':[gene1], 'y':[gene2]})
df_genes=left_gene(df_genes,'x','y',bedfile=args.bedfile,left_col='left')


for i,row in df_genes.iterrows():
    # pedantic gene left/right assignment
    g1=row['x']
    g2=row['y']
    gl=row['left']
    gr= (set([g1,g2]) - set([gl])).pop()
    print(f"# left gene: {gl}, right gene: {gr}")
    # get gene coords
    gl_bed_row=df_bed[df_bed['gene']==gl]
    gr_bed_row=df_bed[df_bed['gene']==gr]
    gl_chrom=gl_bed_row['chrom'].values[0]
    gl_start=gl_bed_row['start'].values[0]
    gl_end=gl_bed_row['end'].values[0]
    gr_chrom=gr_bed_row['chrom'].values[0]
    gr_start=gr_bed_row['start'].values[0]
    gr_end=gr_bed_row['end'].values[0]
    # add left/right gene bed info to stix summary
    df_f['left']=gl
    df_f['right']=gr
    df_f['left_gene_chrom']=gl_chrom
    df_f['left_gene_start']=gl_start
    df_f['left_gene_end']=gl_end
    df_f['right_gene_chrom']=gr_chrom
    df_f['right_gene_start']=gr_start
    df_f['right_gene_end']=gr_end

# reorder columns to have gene info first
leading=df_f.columns.tolist()[31:]
trailing=df_f.columns.tolist()[:31]
new_col_order=leading+trailing
df_f=df_f[new_col_order]
df_f.to_csv(args.output, sep='\t', index=False)
print(f"# wrote output to {args.output}")
print("# completed in %.2f seconds"%(time.time()-t_0))