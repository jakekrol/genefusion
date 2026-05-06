# see stix2fusion for preprocessing steps
import gzip
import shlex
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from polymerization.io import *
import os
import pandas as pd
import numpy as np
import pysam

def validate_bgzip():
    bgzip = shutil.which('bgzip')
    if not bgzip:
        raise FileNotFoundError("bgzip command not found in PATH")
    return bgzip

def run_giggle(argstring, left_gene, outfile, timeout = 60 * 60 * 2, bgzip=False):
    '''
    run giggle with given arguments and return the output as a pandas dataframe
    argstring: string of arguments to pass to giggle, should include {left_gene} and as placeholders for the gene names
    left_gene: name of left gene
    outfile: path to output file to write giggle results to
    timeout: time limit for giggle to run in seconds, default is 2 hours
    returns a pandas dataframe of giggle results
    '''
    giggle = shutil.which('giggle')
    if not giggle:
        raise FileNotFoundError("giggle command not found in PATH")
    cmd = [giggle] + shlex.split(argstring)
    output_path = outfile if not bgzip or outfile.endswith('.gz') else f"{outfile}.gz"
    if bgzip:
        bgzip_cmd = validate_bgzip()
        try:
            with open(output_path, 'wb') as out_handle:
                with subprocess.Popen([bgzip_cmd, '-c'], stdin=subprocess.PIPE, stdout=out_handle, stderr=subprocess.PIPE) as bgzip_proc:
                    header = f"#gene_left={left_gene}\n".encode()
                    bgzip_proc.stdin.write(header)
                    bgzip_proc.stdin.flush()

                    print(f"Running command: {' '.join(cmd)}")
                    giggle_proc = subprocess.Popen(cmd, stdout=bgzip_proc.stdin, stderr=subprocess.PIPE)
                    _, giggle_stderr = giggle_proc.communicate(timeout=timeout)
                    bgzip_proc.stdin.close()

                    bgzip_stderr = bgzip_proc.communicate()[1]

                    if giggle_proc.returncode != 0:
                        print(f"Giggle command failed with error code {giggle_proc.returncode}")
                        if giggle_stderr:
                            print(f"Stderr: {giggle_stderr.decode()}")
                        raise subprocess.CalledProcessError(giggle_proc.returncode, cmd, stderr=giggle_stderr)
                    if bgzip_proc.returncode != 0:
                        print(f"bgzip command failed with error code {bgzip_proc.returncode}")
                        if bgzip_stderr:
                            print(f"Stderr: {bgzip_stderr.decode()}")
                        raise subprocess.CalledProcessError(bgzip_proc.returncode, [bgzip_cmd, '-c'], stderr=bgzip_stderr)
        except subprocess.TimeoutExpired as e:
            print(f"Giggle command timed out after {timeout} seconds")
            raise e
    else:
        # run the command using subprocess.run with a timeout
        try:
            with open(output_path, 'w') as f:
                # write genes in the header
                f.write(f"#gene_left={left_gene}\n")
                f.flush()
                print(f"Running command: {' '.join(cmd)}")
                subprocess.run(cmd, text=True, check=True, stdout=f, stderr=subprocess.PIPE, timeout=timeout)
        except subprocess.CalledProcessError as e:
            print(f"Giggle command failed with error code {e.returncode}")
            print(f"Stderr: {e.stderr.decode()}")
            raise e
        except subprocess.TimeoutExpired as e:
            print(f"Giggle command timed out after {timeout} seconds")
            raise e
    return output_path

def merge_fusion_set_bed2giggle(
    df_merged,
    df_shard,
    outdir,
    outfile_prefix='',
    outfile_suffix='.giggle',
    max_workers=4,
    gene_delim='--',
    timeout=60 * 60 * 2,
    bgzip=True,
    outfile_column_giggle = 'outfile_giggle'
):
    # similar to stix2fusion.merge_fusion_set_bed2stix() but for giggle instead of stix

    # only left-genes are searched
    # significant reduction in complexity
    genes = set(df_merged['gene_left'])
    max_workers=min(max_workers, len(genes), os.cpu_count())
    print(f"# giggle searching for {len(genes)} genes using {max_workers} workers")

    df_merged[outfile_column_giggle] = '' # initialize outfile column
    for cat,group in df_shard.groupby('category'):
        # each category get its own directory
        outdir_cat = os.path.join(outdir, cat)
        os.makedirs(outdir_cat, exist_ok=True)
        giggle_index = group['giggle_index'].iloc[0]
        # parallel giggle queries
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            # for each gene
            for gene in genes:
                row = df_merged[df_merged['gene_left'] == gene]
                # query data
                gene_left = row['gene_left'].iloc[0]

                chromosome_left = row['chromosome_left'].iloc[0]
                start_left = row['start_left'].iloc[0]
                end_left = row['end_left'].iloc[0]
                region_left = f"{chromosome_left}:{start_left}-{end_left}"

                # construct giggle command arguments
                argstring = f"search -i {giggle_index} -r {region_left} -v"
                outfile = os.path.join(
                    outdir_cat, f"{outfile_prefix}{gene_left}{outfile_suffix}"

                )
                if bgzip and not outfile.endswith('.gz'):
                    outfile += '.gz'
                df_merged.loc[df_merged['gene_left'] == gene, outfile_column_giggle] = outfile
                # run giggle
                ex.submit(run_giggle, argstring, gene_left, outfile, timeout, bgzip)
    # includes giggle outfile column
    df_giggle = df_merged.copy()
    return df_giggle

def swap_intervals(path_gigglefile, outfile, bgzip=False):
    # input columns:
    # 1. chrom_left
    # 2. start_left
    # 3. end_left
    # 4. strand_left
    # 5. chrom_right
    # 6. start_right
    # 7. end_right
    # 8. strand_right
    # 9. evidence type
    # 10. sample id
    # output columns:
    # 1. chrom_right
    # 2. start_right
    # 3. end_right
    # 4. strand_right
    # 5. chrom_left
    # 6. start_left
    # 7. end_left
    # 8. strand_left
    # 9. evidence type
    # 10. sample id
    if bgzip:
        mode = 'rt'
        with gzip.open(path_gigglefile, mode) as f_in:
            with pysam.BGZFile(outfile, 'wb') as f_out:
                for line in f_in:
                    if line.startswith('#'):
                        f_out.write(line.encode())
                    else:
                        fields = line.strip().split('\t')
                        f_out.write('\t'.join([
                            fields[4], fields[5], fields[6], fields[7],
                            fields[0], fields[1], fields[2], fields[3],
                            fields[8], fields[9]
                        ]).encode() + b'\n')
    else:
        mode = 'r'
        with open(path_gigglefile, mode) as f_in:
            with open(outfile, 'w') as f_out:
                for line in f_in:
                    if line.startswith('#'):
                        f_out.write(line)
                    else:
                        fields = line.strip().split('\t')
                        f_out.write('\t'.join([
                            fields[4], fields[5], fields[6], fields[7],
                            fields[0], fields[1], fields[2], fields[3],
                            fields[8], fields[9]
                        ]) + '\n')

def giggle2swap(df_giggle, bgzip=False, outfile_column_giggle='outfile_giggle', outfile_column_swap='outfile_swap'):
    breakpoint()
    df_swap = df_giggle.copy()
    # swap left and right intervals in giggle output
    for idx, row in df_giggle.iterrows():
        infile = row[outfile_column_giggle]
        outfile = infile.replace('.giggle', '.giggle.swap')
        swap_intervals(infile, outfile, bgzip)
        df_swap.at[idx, outfile_column_swap] = outfile
    return df_swap

        
