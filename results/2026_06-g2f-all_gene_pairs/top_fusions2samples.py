#!/usr/bin/env python3

import argparse
from polymerization.io import *
from polymerization.analysis import *
import pandas as pd
import os

parser = argparse.ArgumentParser(description="Get samples of top fusions")
parser.add_argument("--input", "-i", required=True, help="Input top fusion tsv")
parser.add_argument("--tissue", "-t", required=True, help="Tissue name")
parser.add_argument("--g2f_outdir", "-g", default="./g2f_out", help="Output directory of g2f results")
parser.add_argument("--output", "-o", required=True, help="Output tsv of top fusions with samples")
parser.add_argument('--evidence_file_suffix', default='.giggle.clean.swap.intersect.bed.gz', help='Suffix for evidence file')
parser.add_argument('--bed', default='/data/jake/genefusion/results/2025_04-gene_bedfile_cln/grch37.genes.bed')
parser.add_argument('--id2name_map', default='/data/jake/genefusion/results/2025_05-pcawg_filename2id/filename2fileid.tsv')
parser.add_argument('--id2specimen_map', default='/data/jake/genefusion/results/2025_05-pcawg_fileid2sample_type/fileid2sampletype.tsv')
args = parser.parse_args()

def main():
	df_name2id = pd.read_csv(args.id2name_map, sep='\t')
	df_specimen = pd.read_csv(args.id2specimen_map, sep='\t')
	df_fusion = pd.read_csv(args.input, sep="\t")
	df_bed = read_bed(args.bed, gene_col_idx=3)
	dir_tissue_g2f = os.path.join(args.g2f_outdir, f"{args.tissue}_tumor_rna")
	n = len(df_fusion)
	# df_fusion['samples'] = ''
	# df_fusion['left_chromosome'] = ''
	# df_fusion['right_chromosome'] = ''
	# df_fusion['left_start'] = -1
	# df_fusion['left_end'] = -1
	# df_fusion['right_start'] = -1
	# df_fusion['right_end'] = -1
	with open(args.output, 'w') as f_out:
		f_out.write("fusion\tgene_left\tgene_right\tscore_uniform\tannotation\tsample\tfilename\tspecimen\tleft_chromosome\tright_chromosome\tleft_breakpoint\tright_breakpoint\n")
		for i, row in df_fusion.iterrows():
			print("# processing {}/{}: {}-{}".format(i+1, n, row['gene_left'], row['gene_right']))
			gene_left = row['gene_left']
			gene_right = row['gene_right']
			left_chromosome = df_bed.loc[df_bed['gene_name'] == gene_left, 'chromosome'].values[0]
			right_chromosome = df_bed.loc[df_bed['gene_name'] == gene_right, 'chromosome'].values[0]
			evidence_file = os.path.join(dir_tissue_g2f, f"{gene_left}{args.evidence_file_suffix}")
			_, df_intersect = read_g2f_intersect(evidence_file, sep='\t', header=None, bgzip=True)
			df_breakpoints = intersect2breakpoints(df_intersect, gene_right, group_by_sample=True)
			df_breakpoints['sample'] = df_breakpoints['sample'].apply(lambda x: os.path.basename(x).split('.')[0])
			for j, rowbp in df_breakpoints.iterrows():
				sample = rowbp['sample']
				# map sample to fileid
				mask = df_name2id['File_ID'] == sample
				filename = df_name2id[mask]['File_Name'].values[0] if mask.any() else ''
				# map fileid to specimen
				mask = df_specimen['File_ID'] == sample
				specimen = df_specimen[mask]['Specimen_Type'].values[0] if mask.any() else ''
				f_out.write("{}--{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
					gene_left, gene_right, row['gene_left'], row['gene_right'], row['score_uniform'], row['annotation'],
					sample, filename, specimen, left_chromosome, right_chromosome,
					rowbp['left_breakpoint'], rowbp['right_breakpoint']
				))


if __name__ == "__main__":
	main()