#!/usr/bin/env python3
import argparse
import ast
import glob
import os
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("--exons", default= "../2026_09-grch37-gtf2exon_intervals/exon_intervals.tsv.gz")
parser.add_argument("--indir",default="top_fusions")
parser.add_argument("--infile_suffix",default = ".annotate.repeat_filt.breakpoint.tsv")
parser.add_argument("--outfile_suffix", default = ".annotate.repeat_filt.breakpoint.rank.tsv")
parser.add_argument("--outdir", default="top_fusions")
parser.add_argument("--default_score",default=2**32)
args = parser.parse_args()

def validate_args():
    assert os.path.exists(args.exons)
    assert os.path.isdir(args.indir)
    assert os.path.isdir(args.outdir)

def score_breakpoint(gene, df_exon, bp,default_score):
    df_sub = df_exon.loc[[gene], :] # wrapping gene in a list ensure result is df, not series which occured sometimes
    score = default_score
    transcript = None
    try:
        for i,row in df_sub.iterrows():
            transcript_name = row['transcript_name']
            exon_intervals = row['exon_intervals']
            exon_intervals = ast.literal_eval(exon_intervals) # string to list of tuples
            for exon_interval in exon_intervals:
                start=int(exon_interval[0])
                end=int(exon_interval[1])
                new_score = min(
                    abs(bp-start),
                    abs(bp-end)
                )
                if new_score < score:
                    score = new_score
                    transcript = transcript_name
    except AttributeError:
        breakpoint()
    return transcript,score

def main():
    df_exon = pd.read_csv(args.exons,sep="\t")
    df_exon.set_index("gene",inplace=True)
    files = glob.glob(f"{args.indir}/*{args.infile_suffix}")
    for f in files:
        fbase = os.path.basename(f)
        fout = fbase.replace(args.infile_suffix,args.outfile_suffix)
        fout = os.path.join(
            args.outdir,
            fout
        )
        df = pd.read_csv(f,sep="\t")
        df['breakpoint_rank'] = args.default_score
        df['breakpoint_rank'] = df['breakpoint_rank'].astype(float)
        df['transcript_left'] = pd.NA
        df['transcript_right'] = pd.NA
        for i, row in df.iterrows():
            gene_left = row['gene_left']
            gene_right = row['gene_right']
            bp_left = row['breakpoint_left']
            bp_right = row['breakpoint_right']
            if (pd.isna(bp_left) or pd.isna(bp_right)):
                continue
            bp_left = int(bp_left.split(":")[1])
            bp_right = int(bp_right.split(":")[1])
            transcript_name_left,left_bp_score = score_breakpoint(gene_left,df_exon,bp_left, args.default_score)
            transcript_name_right,right_bp_score = score_breakpoint(gene_right,df_exon,bp_right, args.default_score)
            df.at[i, 'transcript_left'] = transcript_name_left
            df.at[i, 'transcript_right'] = transcript_name_right
            bp_score_avg = (left_bp_score + right_bp_score) / 2
            df.at[i, 'breakpoint_rank'] = bp_score_avg
        df = df.sort_values('breakpoint_rank')
        df.to_csv(fout,sep='\t',index=False)


if __name__ == "__main__":
    main()