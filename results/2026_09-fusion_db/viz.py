#!/usr/bin/env python3
import argparse
import ast
import glob
import os
import pandas as pd
import subprocess
import tempfile
parser = argparse.ArgumentParser()
parser.add_argument("--pcawg_file_id2ega_file_id", default="../2026_10-pcawg2ega/pcawg_fileid2ega_fileid.tsv")
parser.add_argument("--dir_bam",default="bams")
parser.add_argument("--dir_table", default="top_fusions")
parser.add_argument("--infile_suffix", default=".annotate.repeat_filt.breakpoint.rank.samples.tsv")
parser.add_argument("--outdir_plot", default='samplots')
parser.add_argument("--outdir_table", default='top_fusions')
parser.add_argument("--outfile_plot_suffix", default=".samplot.png")
parser.add_argument("--outfile_table_suffix", default=".annotate.repeat_filt.breakpoint.rank.samples.plot.tsv")
parser.add_argument("--transcripts", default="../2026_09-grch37-gtf2exon_intervals/exon_intervals.tsv.gz")
parser.add_argument("--conda_samplot", default="samplot_env",help="name of samplot conda env")
args = parser.parse_args()
os.makedirs(args.outdir_plot, exist_ok=True)
os.makedirs(args.outdir_table,exist_ok=True)
assert os.path.exists(args.pcawg_file_id2ega_file_id)
assert os.path.exists(args.transcripts)
assert os.path.isdir(args.dir_bam)
assert os.path.isdir(args.dir_table)

# use transcript as annotation track for samplot
def make_exon_bed(transcript_left, chromosome_left,transcript_right, chromosome_right, df_transcripts, outfile):
    try:
        exons_left = df_transcripts.loc[transcript_left,"exon_intervals"]
    except KeyError as e:
        print(e)
        return None
    try:
        exons_right = df_transcripts.loc[transcript_right,"exon_intervals"]
    except KeyError as e:
        print(e)
        return None
    exons_left = ast.literal_eval(exons_left)
    exons_right = ast.literal_eval(exons_right)
    with open(outfile, "w") as f:
        for interval in exons_left:
            f.write(f"{str(chromosome_left)}\t{interval[0]}\t{interval[1]}\n")
        for interval in exons_right:
            f.write(f"{str(chromosome_right)}\t{interval[0]}\t{interval[1]}\n")
    cmd = f"bgzip -f {outfile}"
    print(f"# running cmd: {cmd}")
    subprocess.run(cmd,shell=True,check=True)
    cmd = f"tabix -f -p bed {outfile}.gz"
    print(f"# running cmd: {cmd}")
    subprocess.run(cmd,shell=True, check=True)
    return f"{outfile}.gz"

def samplot(
    bams,
    name,
    chromosome1,
    start1,
    end1,
    outfile,
    annotation_bed,
    svtype="BND",
    chromosome2=None,
    start2=None,
    end2=None,
    conda_env=args.conda_samplot
):
    cmd = f"samplot plot -b {bams} -n {name} -c {chromosome1} "
    cmd += f"-s {start1} -e {end1} -t {svtype} -A {annotation_bed} "
    cmd += f"-o {outfile} --zoom {10000}"
    if all([chromosome2, start2, end2]):
        cmd += f" -c {chromosome2} -s {start2} -e {end2}"
    if conda_env:
        cmd = f"source ~/.bashrc; conda activate {conda_env}; " + cmd
    print(f"# running cmd: {cmd}")
    result=subprocess.run(cmd,shell=True,check=True)
    return result.returncode


def main():
    df_pcawg2ega = pd.read_csv(args.pcawg_file_id2ega_file_id,sep="\t")
    df_pcawg2ega.set_index("File_ID",inplace=True)
    df_transcripts = pd.read_csv(args.transcripts,sep="\t")
    df_transcripts.set_index("transcript_name",inplace=True)
    tables = glob.glob(f"{args.dir_table}/*{args.infile_suffix}")
    bams = glob.glob(f"{args.dir_bam}/*/*.bam")
    for f in tables:
        df = pd.read_csv(f,sep="\t")
        for i,row in df.iterrows():
            gene_left = row['gene_left']
            gene_right = row['gene_right']
            chromosome_left = row['breakpoint_left'].split(":")[0]
            chromosome_right = row['breakpoint_right'].split(":")[0]
            start = row['breakpoint_left'].split(":")[1]
            end = row['breakpoint_right'].split(":")[1]
            transcript_left=row['transcript_left']
            transcript_right=row['transcript_right']
            samples = row['top_samples'].split(",")
            plot_bams = []
            for sample in samples:
                ega_file_id = df_pcawg2ega.loc[sample, "ega_file_id"]
                for bam in bams:
                    if ega_file_id in bam:
                        # a checksum file means download verification passed
                        if os.path.exists(bam + '.md5'):
                            plot_bams.append(bam)
            if plot_bams:
                tmp_bed = f"./tmp/{transcript_left}--{transcript_right}.bed"
                annotate_bed =  make_exon_bed(
                    transcript_left,
                    chromosome_left,
                    transcript_right,
                    chromosome_right,
                    df_transcripts,
                    tmp_bed
                )
                if annotate_bed:
                    bam_str=','.join(plot_bams)
                    name = f"{gene_left}--{gene_right}"
                    outfile_samplot=os.path.join(
                        args.outdir_plot,
                        f"{gene_left}--{gene_right}{args.outfile_plot_suffix}"
                    )
                    if chromosome_left == chromosome_right:
                        samplot(
                            bam_str,
                            name,
                            chromosome_left,
                            start,
                            end,
                            outfile_samplot,
                            annotate_bed,
                            svtype="BND",
                            chromosome2=None,
                            start2=None,
                            end2=None
                        )
                    else:
                        samplot(
                            bam_str,
                            name,
                            chromosome_left,
                            start,
                            end,
                            outfile_samplot,
                            annotate_bed,
                            svtype="BND",
                            chromosome2=chromosome_right,
                            start2=start,
                            end2=end,
                        )

                            



if __name__ == "__main__":
    main()