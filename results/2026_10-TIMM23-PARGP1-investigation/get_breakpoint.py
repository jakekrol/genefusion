#!/usr/bin/env python3
import argparse
import ast
import os
import pandas as pd
import subprocess
from polymerization.analysis import intersect_file2breakpoints
parser = argparse.ArgumentParser()
parser.add_argument("--samples", default="samples.tsv")
parser.add_argument("--exons", default= "../2026_09-grch37-gtf2exon_intervals/exon_intervals.tsv.gz")
parser.add_argument("--default_score",default=2**32)
parser.add_argument(
    "--intersect_thousg_rna",
    default="../2026_06-g2f-all_gene_pairs/g2f_out/mage_short_read_1000g_rna/PARGP1.giggle.clean.swap.intersect.bed.gz"
)
parser.add_argument(
    "--intersect_pcawg_blood_tumor",
    default="../2026_06-g2f-all_gene_pairs/g2f_out/blood_tumor_rna/PARGP1.giggle.clean.swap.intersect.bed.gz"
)
# parser.add_argument("--dir_thousg_rna",default="../2026_07-1kg-rna/star_fusion-parallel",help="thousg rna bams live in subdirs of here")
parser.add_argument("--output", default="breakpoints.tsv")
args = parser.parse_args()

GENE_LEFT='PARGP1'
GENE_RIGHT='TIMM23'

def get_subset_intersect(file_intersect: str, samples: list):
    df = intersect_file2breakpoints(file_intersect, GENE_RIGHT, False, True)
    mask = df['sample'].isin(samples)
    df = df[mask].reset_index(drop=True)
    return df
    

### borrowed from /data/jake/genefusion/results/2026_09-fusion_db/estimate_breakpoint.py
def estimate_breakpoint(df_intersect: pd.DataFrame):
    left_chrom = df_intersect.loc[0,'chromosome_left']
    right_chrom = df_intersect.loc[0,'chromosome_right']
    bp_left = int(df_intersect['start_left'].mean())
    bp_right = int(df_intersect['end_right'].mean())
    bp_left = f"{left_chrom}:{bp_left}"
    bp_right = f'{right_chrom}:{bp_right}'
    return bp_left, bp_right

### borrowed from /data/jake/genefusion/results/2026_09-fusion_db/rank_breakpoints.py
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
    mask = df_exon['gene'].isin([GENE_LEFT,GENE_RIGHT])
    df_exon = df_exon[mask].reset_index(drop=True)
    df_exon.set_index('gene',inplace=True)
    df_samples = pd.read_csv(args.samples,sep="\t")
    ### thousg rna
    samples_thousg_rna = df_samples['sample'].iloc[0:3].tolist()
    df_intersect_subset = get_subset_intersect(args.intersect_thousg_rna, samples_thousg_rna)
    bp_left, bp_right = estimate_breakpoint(df_intersect_subset)
    df_samples['breakpoint_left'] = ''
    df_samples['transcript_left'] = ''
    df_samples['breakpoint_right'] = ''
    df_samples['transcript_right'] = ''
    df_samples.loc[0:2,'breakpoint_left'] = bp_left
    df_samples.loc[0:2,'breakpoint_right'] = bp_right
    transcript_left,_ = score_breakpoint(GENE_LEFT, df_exon, int(bp_left.split(':')[1]),args.default_score)
    transcript_right,_ = score_breakpoint(GENE_RIGHT, df_exon, int(bp_right.split(':')[1]),args.default_score)
    df_samples.loc[0:2,'transcript_left'] = transcript_left
    df_samples.loc[0:2,'transcript_right'] = transcript_right
    ### pcawg blood tumor rna
    samples_pcawg_blood_tumor = df_samples['sample'].iloc[3:].tolist()
    df_intersect_subset = get_subset_intersect(args.intersect_pcawg_blood_tumor, samples_pcawg_blood_tumor)
    bp_left, bp_right = estimate_breakpoint(df_intersect_subset)
    df_samples.loc[3:,'breakpoint_left'] = bp_left
    df_samples.loc[3:,'breakpoint_right'] = bp_right
    transcript_left,_ = score_breakpoint(GENE_LEFT, df_exon, int(bp_left.split(':')[1]),args.default_score)
    transcript_right,_ = score_breakpoint(GENE_RIGHT, df_exon, int(bp_right.split(':')[1]),args.default_score)
    df_samples.loc[3:,'transcript_left'] = transcript_left
    df_samples.loc[3:,'transcript_right'] = transcript_right
    df_samples.to_csv(args.output,sep="\t",index=False)

if __name__ == "__main__":
    main()

# ### borrowed from /data/jake/genefusion/results/2026_09-fusion_db/viz.py
# # use transcript as annotation track for samplot
# def make_exon_bed(transcript_left, chromosome_left,transcript_right, chromosome_right, df_transcripts, outfile):
#     try:
#         exons_left = df_transcripts.loc[transcript_left,"exon_intervals"]
#     except KeyError as e:
#         print(e)
#         return None
#     try:
#         exons_right = df_transcripts.loc[transcript_right,"exon_intervals"]
#     except KeyError as e:
#         print(e)
#         return None
#     exons_left = ast.literal_eval(exons_left)
#     exons_right = ast.literal_eval(exons_right)
#     with open(outfile, "w") as f:
#         for interval in exons_left:
#             f.write(f"{str(chromosome_left)}\t{interval[0]}\t{interval[1]}\n")
#         for interval in exons_right:
#             f.write(f"{str(chromosome_right)}\t{interval[0]}\t{interval[1]}\n")
#     # sort
#     cmd = f"bedtools sort -i {outfile} > {outfile.replace(".bed",".sort.bed")}"
#     print(f"# running cmd: {cmd}")
#     subprocess.run(cmd,shell=True,check=True)
#     outfile = f"{outfile.replace(".bed",".sort.bed")}"
#     # compress
#     cmd = f"bgzip -f {outfile}"
#     print(f"# running cmd: {cmd}")
#     subprocess.run(cmd,shell=True,check=True)
#     outfile = f"{outfile}.gz"
#     # index
#     cmd = f"tabix -f -p bed {outfile}"
#     print(f"# running cmd: {cmd}")
#     subprocess.run(cmd,shell=True, check=True)
#     return f"{outfile}"

# def samplot(
#     bams,
#     name,
#     chromosome1,
#     start1,
#     end1,
#     outfile,
#     annotation_bed,
#     svtype="BND",
#     chromosome2=None,
#     start2=None,
#     end2=None,
#     conda_env=args.conda_samplot
# ):
#     cmd = f"samplot plot -b {bams} -n {name} -c {chromosome1} "
#     cmd += f"-s {start1} -e {end1} -t {svtype} -A {annotation_bed} "
#     cmd += f"-o {outfile} --zoom {10000}"
#     if all([chromosome2, start2, end2]):
#         cmd += f" -c {chromosome2} -s {start2} -e {end2}"
#     if conda_env:
#         cmd = f"source ~/.bashrc; conda activate {conda_env}; " + cmd
#     print(f"# running cmd: {cmd}")
#     result=subprocess.run(cmd,shell=True,check=True)
#     return result.returncode