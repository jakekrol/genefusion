#!/usr/bin/env python3

# use samplot_env to have samtools view from aws work
import argparse
import subprocess
import os
import pandas as pd

parser = argparse.ArgumentParser(description="Download region BAMs from aws")
parser.add_argument("--input", "-i", required=True, help="Input tsv of top fusions with samples")
parser.add_argument("--output_dir", "-o", required=True, help="Output directory for downloaded BAMs")
parser.add_argument("--bucket", "-b", required=True, help="S3 bucket name")
parser.add_argument("--padding", "-p", type=int, default=1e6, help="Padding around breakpoints for BAM download")
parser.add_argument("--tissue", "-t", required=True, help="Tissue name")
parser.add_argument("--log", "-l", default="download_region_bams.log", help="Log file")
parser.add_argument("--out_tbl", "-ot", required=True)
args = parser.parse_args()


def main():
	dir_tissue=os.path.join(args.output_dir, args.tissue)
	os.makedirs(dir_tissue, exist_ok=True)
	df = pd.read_csv(args.input, sep='\t')
	df['region_bam'] = ''
	fail=0
	success=0
	with open(args.log, 'w') as log_f:
		for i, row in df.iterrows():
			fusion = row['fusion']
			score = row['score_uniform']
			sample = row['sample']
			filename = row['filename']
			specimen = row['specimen']
			left_chromosome = row['left_chromosome']
			right_chromosome = row['right_chromosome']
			left_breakpoint = row['left_breakpoint']
			right_breakpoint = row['right_breakpoint']
			left_start = int(max(0, left_breakpoint - args.padding))
			left_end = int(left_breakpoint + args.padding)
			right_start = int(max(0, right_breakpoint - args.padding))
			right_end = int(right_breakpoint + args.padding)
			region1 = f"{left_chromosome}:{left_start}-{left_end}"
			region2 = f"{right_chromosome}:{right_start}-{right_end}"
			outfile=os.path.join(dir_tissue, f"{fusion}.{sample}.{specimen}.bam")
			df.at[i, 'region_bam'] = outfile
			try:
				# download index
				cmd = [
					"aws", "s3", "cp",
					f"s3://{args.bucket}/rna_cancerdata/{args.tissue}/{filename}.bai",
					f"{outfile}.bai"
				]
				print("# running command: {}".format(' '.join(cmd)))
				result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				log_f.write(f"# command: {' '.join(cmd)}\n")
				log_f.write(f"# return code: {result.returncode}\n")
				log_f.write(f"# stdout: {result.stdout.decode('utf-8')}\n")
				log_f.write(f"# stderr: {result.stderr.decode('utf-8')}\n")
				cmd = [
					"samtools", "view", "-b", "-o", outfile,
					f"s3://{args.bucket}/rna_cancerdata/{args.tissue}/{filename}",
					f"{region1}", f"{region2}"
				]
				print("# running command: {}".format(' '.join(cmd)))
				result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				# log the result
				retcode = result.returncode
				stdout = result.stdout.decode('utf-8')
				stderr = result.stderr.decode('utf-8')
				log_f.write(f"# command: {' '.join(cmd)}\n")
				log_f.write(f"# return code: {retcode}\n")
				log_f.write(f"# stdout: {stdout}\n")
				log_f.write(f"# stderr: {stderr}\n")
				success += 1
			except subprocess.CalledProcessError as e:
				log_f.write(f"# command: {' '.join(e.cmd)}\n")
				log_f.write(f"# return code: {e.returncode}\n")
				log_f.write(f"# stdout: {e.stdout.decode('utf-8') if e.stdout else ''}\n")
				log_f.write(f"# stderr: {e.stderr.decode('utf-8') if e.stderr else ''}\n")
				log_f.write(f"# error downloading BAM for fusion {fusion}, sample {sample}, specimen {specimen}\n")
				fail += 1
		log_f.write(f"# total successful downloads: {success}\n")
		log_f.write(f"# total failed downloads: {fail}\n")


if __name__ == "__main__":
    main()