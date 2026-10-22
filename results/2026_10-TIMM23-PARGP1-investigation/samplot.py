#!/usr/bin/env python3

import argparse
import ast
import os
import pandas as pd
import subprocess
import tempfile
# from polymerization.plot import samplot

parser=argparse.ArgumentParser()
parser.add_argument("--input",default="breakpoints.tsv")
parser.add_argument("--out_normal",default="PARGP1_TIMM23-thousg_rna.png")
parser.add_argument("--out_tumor",default="PARGP1_TIMM23-pcawg_blood_tumor_rna.png")
parser.add_argument("--transcripts", default="../2026_09-grch37-gtf2exon_intervals/exon_intervals.tsv.gz")
parser.add_argument("--conda_samplot", default="samplot_env",help="name of samplot conda env")
parser.add_argument("--dir_thousg_rna",default="../2026_07-1kg-rna/star_fusion-parallel",help="thousg rna bams live in subdirs of here")
parser.add_argument("--dir_bam_out",default="bam")
parser.add_argument("--threads",default=8)
parser.add_argument("--annotation",default="exons.sort.bed.gz")
parser.add_argument("--zoom", default=70_000)
parser.add_argument("--pcawg2ega",default="../2026_10-pcawg2ega/pcawg_fileid2ega_fileid.tsv")
args = parser.parse_args()
args.threads = int(args.threads)

BAM_STAR_FUSION='Aligned.out.bam'
GENE_LEFT='PARGP1'
GENE_RIGHT='TIMM23'

def sample2bamthousg(sample: str, dir_thousg_rna: str =args.dir_thousg_rna):
    sample = os.path.basename(sample)
    srr_id = sample.split('.')[1]
    bam = os.path.join(
        dir_thousg_rna,
        srr_id,
        BAM_STAR_FUSION
    )
    assert os.path.exists(bam)
    return bam

def fid2bampcawg(fid: str):
    # we assume the bams are in subdirs named by fid; standard naming for ega downloads
    files = os.listdir(fid)
    for f in files:
        if f.endswith('.bam'):
            return os.path.join(fid,f)

    

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
    # sort
    cmd = f"bedtools sort -i {outfile} > {outfile.replace(".bed",".sort.bed")}"
    print(f"# running cmd: {cmd}")
    subprocess.run(cmd,shell=True,check=True)
    outfile = f"{outfile.replace(".bed",".sort.bed")}"
    # compress
    cmd = f"bgzip -f {outfile}"
    print(f"# running cmd: {cmd}")
    subprocess.run(cmd,shell=True,check=True)
    outfile = f"{outfile}.gz"
    # index
    cmd = f"tabix -f -p bed {outfile}"
    print(f"# running cmd: {cmd}")
    subprocess.run(cmd,shell=True, check=True)
    return f"{outfile}"

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
    conda_env=args.conda_samplot,
    zoom=args.zoom
):
    cmd = f"samplot plot -b {bams} -n {name} -c {chromosome1} "
    cmd += f"-s {start1} -e {end1} -t {svtype} -A {annotation_bed} "
    cmd += f"-o {outfile} --zoom {zoom}"
    if all([chromosome2, start2, end2]):
        cmd += f" -c {chromosome2} -s {start2} -e {end2}"
    if conda_env:
        cmd = f"source ~/.bashrc; conda activate {conda_env}; " + cmd
    print(f"# running cmd: {cmd}")
    result=subprocess.run(cmd,shell=True,check=True)
    return result.returncode

def prep_bam(bam_in: str, bam_out: str):
    cmd = f"samtools sort -@ {args.threads} -o {bam_out} {bam_in}"
    print(f"# running cmd: {cmd}")
    result=subprocess.run(cmd,shell=True,check=True)
    cmd = f"samtools index -@ {args.threads} {bam_out}"
    print(f"# running cmd: {cmd}")
    result=subprocess.run(cmd,shell=True,check=True)
    return result.returncode

def main():
    df = pd.read_csv(args.input,sep="\t")
    df['bam'] = ''
    df.loc[0:2,'bam'] = df.loc[0:2,'sample'].apply(sample2bamthousg)
    pcawg2ega = pd.read_csv(args.pcawg2ega,sep="\t").set_index("File_ID").squeeze()
    df['ega']=''
    df['pcawg_fid'] = ''
    df.loc[3:,'ega'] = df.loc[3:,'sample'].apply(lambda x: os.path.basename(x).split('.')[0])
    df.loc[3:,'pcawg_fid'] = df.loc[3:,'ega'].map(pcawg2ega)
    df.loc[3:,'bam'] = df.loc[3:,'pcawg_fid'].map(fid2bampcawg)
    df_transcripts = pd.read_csv(args.transcripts,sep="\t")
    df_transcripts = df_transcripts.set_index("transcript_name")
    ### thousg rna
    # sort and index bams
    for i,row in (df.loc[0:2,:]).iterrows():
        sample=os.path.basename(row['sample']).replace('.bed.gz', '.bam')
        bam_in = row['bam']
        bam_out = os.path.join(args.dir_bam_out, sample)
        if not os.path.exists(bam_out):
            prep_bam(bam_in,bam_out)
        df.loc[i,'bam'] = bam_out
    name='1KG-MAGE-PARGP1_TIMM23'
    bam_str = ' '.join(df.loc[0:2,'bam'].tolist())
    transcript_left = df.loc[0,'transcript_left']
    transcript_right = df.loc[0,'transcript_right']
    chromosome_left = df.loc[0,'breakpoint_left'].split(':')[0]
    chromosome_right = df.loc[0,'breakpoint_right'].split(':')[0]
    start = df.loc[0,'breakpoint_left'].split(':')[1]
    end = df.loc[0,'breakpoint_right'].split(':')[1]
    # tmp_bed = f"./tmp/{transcript_left}--{transcript_right}.bed"
    # annotation_bed = make_exon_bed(transcript_left, chromosome_left,transcript_right, chromosome_right, df_transcripts, tmp_bed)
    samplot(
        bam_str,
        name,
        chromosome_left,
        start,
        end,
        args.out_normal,
        args.annotation,
        svtype="BND",
        chromosome2=None,
        start2=None,
        end2=None
    )
    ### pcawg blood tumor rna
     # sort and index bams
    for i,row in (df.loc[3:,:]).iterrows():
        fid=os.path.basename(row['sample']).replace('.bed.gz', '.bam')
        bam_in = row['bam']
        bam_out = os.path.join(args.dir_bam_out, fid)
        if not os.path.exists(bam_out):
            prep_bam(bam_in,bam_out)
        df.loc[i,'bam'] = bam_out
    name='PCAWG_BLOOD_TUMOR_RNA-PARGP1_TIMM23'
    bam_str = ' '.join(df.loc[3:,'bam'].tolist())
    transcript_left = df.loc[3,'transcript_left']
    transcript_right = df.loc[3,'transcript_right']
    chromosome_left = df.loc[3,'breakpoint_left'].split(':')[0]
    chromosome_right = df.loc[3,'breakpoint_right'].split(':')[0]
    start = df.loc[3,'breakpoint_left'].split(':')[1]
    end = df.loc[3,'breakpoint_right'].split(':')[1]
    samplot(
        bam_str,
        name,
        chromosome_left,
        start,
        end,
        args.out_tumor,
        args.annotation,
        svtype="BND",
        chromosome2=None,
        start2=None,
        end2=None,
        zoom=args.zoom
    )

if __name__ == "__main__":
    main()