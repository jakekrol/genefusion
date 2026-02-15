#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import pandas as pd
import time

t_0=time.time()
parser=argparse.ArgumentParser(description="Plot all the results")
parser.add_argument('-i',"--input", help="Input directory with all the results", required=True)
parser.add_argument('-o', "--outdir", help="Output directory for the plots", required=True)
parser.add_argument('-b', '--bucket', help='S3 bucket to download the files from', required=True)
# Homo_sapiens.GRCh37.82.gene_only.sort.gff3.gz
parser.add_argument('-t','--transcript_file', help='Transcript file in GTF format', 
                    default=None)
# grch37.genes.promoter_pad.bed.gz.tbi
parser.add_argument('-a','--annotation_file', help='Annotation bed file',
                    default=None)
parser.add_argument('-c', '--cache', help='Use cached results', action='store_true')
args=parser.parse_args()
if not os.environ.get('AWS_SECRET_ACCESS_KEY') or not os.environ.get('AWS_ACCESS_KEY_ID'):
    print('warning: AWS credentials not found')
else:
    print('AWS credentials found')
# make out dir
os.makedirs(args.outdir, exist_ok=True)
# make out subdir for bam regions and indices
os.makedirs(f'{args.outdir}/bams', exist_ok=True)
# make out subdir for plots
os.makedirs(f'{args.outdir}/plots', exist_ok=True)
if args.annotation_file:
    os.makedirs(f'{args.outdir}/annotations', exist_ok=True)

df = pd.read_csv(args.input, sep="\t")
groups = df.groupby(['tissue','svtype','left','right'])
total_groups = groups.ngroups
group_counter = 0

# for each fusion
for g, df_g in groups:
    group_counter += 1
    print(f"Processing group {group_counter} of {total_groups}: {g[0]}_{g[2]}--{g[3]}_{g[1]}")
    subplots=[]
    evidence_across_samples=0
    # for each sample supporting the fusion
    for i,row in df_g.iterrows():
        leftgene=row['left']
        rightgene=row['right']
        chrom=row['left_gene_chrom']
        start=row['left_breakpoint']
        end=row['right_breakpoint']
        svtype=row['svtype']
        rnafile=row['rna_file_name']
        region=row['region_string']
        total_read_evidence=row['tumor_paired_plus_split_samplewise_evidence']
        tissue=row['tissue']
        project=row['Project']
        donor_id=row['ICGC_Donor']
        specimen_type=row['Specimen_Type']
        specimen_type=specimen_type.replace(' ', '_')
        specimen_type=specimen_type.replace('-', '_')
        # replace repeats of _ with single _
        while '__' in specimen_type:
            specimen_type=specimen_type.replace('__', '_')

        awspath = f's3://{args.bucket}/rna_cancerdata/{tissue}/{rnafile}'
        bamfileout=f'{args.outdir}/bams/{tissue}_{leftgene}--{rightgene}_{svtype}_{rnafile}.bam'
        if bamfileout.endswith('.bam.bam'):
            bamfileout = bamfileout[:-4]
        name_string=f'{donor_id}-{rnafile}-{specimen_type}-ev{total_read_evidence}'
        if (args.cache) and (os.path.exists(bamfileout)) and (os.path.exists(bamfileout+'.bai')):
            print(f"Using cached BAM for {tissue}_{leftgene}--{rightgene}_{svtype} from {rnafile}")
            subplots.append((bamfileout, name_string))
            evidence_across_samples += total_read_evidence
            continue
        try:
            # download
            download_cmd = f'samtools view -b {awspath} {region} -o {bamfileout}'
            subprocess.run(download_cmd, shell=True, check=True)
            index_cmd = f'samtools index {bamfileout}'
            subprocess.run(index_cmd, shell=True, check=True)
            print(f"Successfully processed {tissue}_{leftgene}--{rightgene}_{svtype}")
            subplots.append((bamfileout, name_string))
            evidence_across_samples += total_read_evidence
        except subprocess.CalledProcessError:
            print(f"Skipping {tissue}_{leftgene}--{rightgene}_{svtype}: file not found or processing failed")
            continue
    
    # plotting
    if not subplots:
        print(f"No BAM files found for {tissue}_{leftgene}--{rightgene}_{svtype}, skipping plot")
        continue
    else:
        bam_string = ' '.join([b[0] for b in subplots])
        name_string = ' '.join([b[1] for b in subplots])
        outfile_plot = f'{args.outdir}/plots/{tissue}_{leftgene}--{rightgene}_{svtype}_ev{evidence_across_samples}.png'
        cmd_plot = f'samplot plot -b {bam_string} -c {chrom} -s {start} -e {end} '\
                f'-o {outfile_plot} '\
                f' -n {name_string}'
                    # f'-t {svtype} -T {args.transcript_file}' \
        if args.transcript_file:
            cmd_plot += f' -T {args.transcript_file}'
        if args.annotation_file:
            # subset gene bed file to only include the gene pairs in the plot
            annotation_file_subset = f'{args.outdir}/annotations/{leftgene}--{rightgene}.bed'
            cmd_make_annotation_subset = f'grep -wE "{leftgene}|{rightgene}" {args.annotation_file} |'\
                f'/data/jake/bedtools.static.binary sort -i stdin | bgzip -c > {annotation_file_subset}.gz && tabix -f {annotation_file_subset}.gz'
            subprocess.run(cmd_make_annotation_subset, shell=True, check=True)
            cmd_plot += f' -A {annotation_file_subset}.gz'
        print(f"# running plot command: {cmd_plot}")
        subprocess.run(cmd_plot, shell=True, check=True)


t_1=time.time()
print(f"# done in {t_1-t_0:.2f} seconds. plots at {args.outdir}/plots")

    


        
    