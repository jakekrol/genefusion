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
import glob

def validate_bgzip():
    bgzip = shutil.which('bgzip')
    if not bgzip:
        raise FileNotFoundError("bgzip command not found in PATH")
    return bgzip

def run_giggle(argstring, left_gene, outfile, timeout = 60 * 60 * 2, bgzip=False, append=False):
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
    write_header = (not append) or (not os.path.exists(output_path)) or (os.path.getsize(output_path) == 0)
    if bgzip:
        bgzip_cmd = validate_bgzip()
        assert outfile.endswith('.gz'), f"Output file {outfile} must end with .gz when bgzip=True"

        try:
            out_mode = 'ab' if append else 'wb'
            with open(output_path, out_mode) as out_handle:
                with subprocess.Popen([bgzip_cmd, '-c'], stdin=subprocess.PIPE, stdout=out_handle, stderr=subprocess.PIPE) as bgzip_proc:
                    if write_header:
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
            out_mode = 'a' if append else 'w'
            with open(output_path, out_mode) as f:
                # write genes in the header once
                if write_header:
                    f.write(f"#gene_left={left_gene}\n")
                    f.flush()
                print(f"Running command: {' '.join(cmd)}")
                subprocess.run(cmd, text=True, check=True, stdout=f, stderr=subprocess.PIPE, timeout=timeout)
        except subprocess.CalledProcessError as e:
            print(f"Giggle command failed with error code {e.returncode}")
            print(f"Stderr: {e.stderr}")
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
    outfile_column_giggle_prefix = 'outfile_giggle'
):
    # similar to stix2fusion.merge_fusion_set_bed2stix() but for giggle instead of stix

    # only left-genes are searched
    # significant reduction in complexity
    genes = set(df_merged['gene_left'])
    max_workers=min(max_workers, len(genes), os.cpu_count())
    print(f"# giggle searching for {len(genes)} genes using {max_workers} workers")
    for i, (cat, group) in enumerate(df_shard.groupby('category')):
        # make outfile column for each category
        col_name = f"{outfile_column_giggle_prefix}_{cat}"
        df_merged[col_name] = '' # initialize outfile column for this category
        # each category get its own directory
        outdir_cat = os.path.join(outdir, cat)
        os.makedirs(outdir_cat, exist_ok=True)
        # loop over different indices in same group
        for giggle_index in group['giggle_index'].unique():
            # parallel giggle queries
            with ThreadPoolExecutor(max_workers=max_workers) as ex:
                futures = []
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
                    df_merged.loc[df_merged['gene_left'] == gene, col_name] = outfile
                    # run giggle
                    futures.append(ex.submit(run_giggle, argstring, gene_left, outfile, timeout, bgzip, append=True))
                for future in as_completed(futures):
                    future.result()
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

def giggle2swap(
    df_giggle,
    df_shard,
    bgzip=False,
    outfile_column_giggle_prefix='outfile_giggle',
    outfile_column_swap_prefix='outfile_swap',
    max_workers=4
):
    df_swap = df_giggle.copy()
    for cat in df_shard['category'].unique():
        col_giggle = f"{outfile_column_giggle_prefix}_{cat}"
        col_swap = f"{outfile_column_swap_prefix}_{cat}"
        df_swap[col_swap] = '' # initialize swap outfile column for this category
        # parallel swap left and right intervals in giggle output
        workers = min(max_workers, os.cpu_count() or 1, len(df_giggle))
        futures = {}
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for idx, row in df_giggle.iterrows():
                infile = row[f"{col_giggle}"]
                outfile = infile.replace('.giggle', '.giggle.swap')
                fut = ex.submit(swap_intervals, infile, outfile, bgzip)
                futures[fut] = (idx, outfile, col_swap)
            for fut in as_completed(futures):
                idx, outfile, col = futures[fut]
                fut.result()
                df_swap.at[idx, col] = outfile
    return df_swap

def validate_bedtools():
    bedtools = shutil.which('bedtools')
    if not bedtools:
        raise FileNotFoundError("bedtools command not found in PATH")
    return bedtools

def bedtools_intersect(path_input, path_bedfile, outfile, gene_col_idx=3, bgzip=False, bedtools_bin=None):
    # output columns:
    # c1-5 are from bed file
    # c6-15 are from swap file
    # validation
    if bedtools_bin:
        bedtools = bedtools_bin
    else:
        bedtools = validate_bedtools()
    assert os.path.exists(path_input), f"Input file {path_input} does not exist"
    assert os.path.exists(path_bedfile), f"Bed file {path_bedfile} does not exist"
    _ = read_bed(path_bedfile, gene_col_idx=gene_col_idx) # validate bed file format
    if bgzip:
        bgzip_cmd = validate_bgzip()
        assert path_input.endswith('.gz'), f"Input file {path_input} must be gzipped when bgzip=True"
        assert outfile.endswith('.gz'), f"Output file {outfile} must end with .gz when bgzip=True"

        with open(outfile, 'wb') as f_out:
            with subprocess.Popen([bgzip_cmd, '-c'], stdin=subprocess.PIPE, stdout=f_out, stderr=subprocess.PIPE) as bgzip_proc:
                # write input header
                with gzip.open(path_input, 'rt') as f_in:
                    for line in f_in:
                        if line.startswith('#'):
                            bgzip_proc.stdin.write(line.encode())
                        else:
                            break
                bgzip_proc.stdin.flush()

                # run bedtools intersect and pipe output to bgzip
                cmd = f"{bedtools} intersect -a <(grep -v '^#' {path_bedfile}) " \
                      f"-b <(zcat {path_input} | grep -v '^#') -wb -wa"
                bedtools_proc = subprocess.Popen(cmd, shell=True, stdout=bgzip_proc.stdin, stderr=subprocess.PIPE)
                bedtools_stderr = bedtools_proc.communicate()[1]
                if bedtools_proc.returncode != 0:
                    print(f"bedtools command failed with error code {bedtools_proc.returncode}")
                    if bedtools_stderr:
                        print(f"Stderr: {bedtools_stderr.decode()}")
                    raise subprocess.CalledProcessError(bedtools_proc.returncode, cmd, stderr=bedtools_stderr)

                bgzip_proc.stdin.close()
                bgzip_stderr = bgzip_proc.communicate()[1]
                if bgzip_proc.returncode != 0:
                    print(f"bgzip command failed with error code {bgzip_proc.returncode}")
                    if bgzip_stderr:
                        print(f"Stderr: {bgzip_stderr.decode()}")
                    raise subprocess.CalledProcessError(bgzip_proc.returncode, [bgzip_cmd, '-c'], stderr=bgzip_stderr)
    else:
        cmd = bedtools + " "\
            f"intersect -a {path_bedfile} " \
            f"-a <(grep -v '^#' {path_bedfile}) " \
            f"-b <(grep -v '^#' {path_input}) " \
            f"-wb -wa"
        try:
            with open(outfile, 'w') as f_out:
                subprocess.run(cmd, shell=True, text=True, check=True, stdout=f_out, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"bedtools command failed with error code {e.returncode}")
            print(f"Stderr: {e.stderr.decode()}")
            raise e

def swap2intersect(
    df_swap,
    df_shard,
    path_bedfile,
    gene_col_idx=3,
    outfile_column_swap_prefix='outfile_swap',
    outfile_column_intersect_prefix='outfile_intersect',
    bgzip=False,
    bedtools_bin=None,
    max_workers=4
):
    df_intersect = df_swap.copy()
    for cat in df_shard['category'].unique():
        col_swap = f"{outfile_column_swap_prefix}_{cat}"
        col_intersect = f"{outfile_column_intersect_prefix}_{cat}"
        df_intersect[col_intersect] = '' # initialize intersect outfile column for this category
        workers = min(max_workers, os.cpu_count() or 1, len(df_swap))
        futures = {}
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for idx, row in df_swap.iterrows():
                infile = row[col_swap]
                if bgzip:
                    outfile = infile.replace('.giggle.swap.gz', '.giggle.swap.intersect.bed.gz')
                else:
                    outfile = infile.replace('.giggle.swap', '.giggle.swap.intersect.bed')
                fut = ex.submit(bedtools_intersect, infile, path_bedfile, outfile, gene_col_idx, bgzip, bedtools_bin)
                futures[fut] = (idx, outfile, col_intersect)
            for fut in as_completed(futures):
                idx, outfile, col = futures[fut]
                fut.result()
                df_intersect.at[idx, col] = outfile
    return df_intersect

    
def intersect2evidence(
    path_intersect,
    outfile,
    right_gene_col=3,
    sample_column=14,
    bgzip=False,
    burden=False
):
    # get left gene from header
    if bgzip:
        with gzip.open(path_intersect, 'rt') as f:
            for line in f:
                if line.startswith('#gene_left='):
                    gene_left = line.strip().split('=')[1]
                    break
    else:
        with open(path_intersect, 'r') as f:
            for line in f:
                if line.startswith('#gene_left='):
                    gene_left = line.strip().split('=')[1]
                    break
    # handles bgzipped or not automatically
    df = pd.read_csv(path_intersect, sep='\t', header=None, comment='#', usecols=[right_gene_col, sample_column])
    df.columns=['gene_right', 'sample_id']
    # read counts
    right_gene_counts = df.groupby('gene_right').size().reset_index(name='reads')
    # sample counts
    sample_counts = df.groupby('gene_right').agg(samples=('sample_id', 'nunique')).reset_index()
    # merge counts
    df_evidence = pd.merge(right_gene_counts, sample_counts, on='gene_right', how='outer')
    df_evidence['gene_left'] = gene_left
    df_evidence = df_evidence[['gene_left', 'gene_right', 'reads', 'samples']]
    if burden:

        burden_outfile = os.path.join(
            os.path.dirname(outfile),
            gene_left + '.burden.txt'
        )
        # rm selfies
        total_burden = df_evidence[df_evidence['gene_left'] != df_evidence['gene_right']]['reads'].sum()
        with open(burden_outfile, 'w') as f:
            f.write(f"{total_burden}\n")
    # write output
    df_evidence.to_csv(outfile, sep='\t', index=False, compression=None)

def df_intersect2df_evidence(
    df_intersect,
    df_shard,
    right_gene_col=3,
    sample_column=14,
    outfile_column_intersect_prefix='outfile_intersect',
    outfile_column_evidence_prefix='outfile_evidence',
    bgzip=False,
    burden=False,
    max_workers=4
):
    df_evidence = df_intersect.copy()
    for cat in df_shard['category'].unique():
        col_intersect = f"{outfile_column_intersect_prefix}_{cat}"
        col_evidence = f"{outfile_column_evidence_prefix}_{cat}"
        df_evidence[col_evidence] = '' # initialize evidence outfile column for this category
        workers = min(max_workers, os.cpu_count() or 1, len(df_intersect))
        futures = {}
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for idx, row in df_intersect.iterrows():
                infile = row[col_intersect]
                if bgzip:
                    outfile = infile.replace('.giggle.swap.intersect.bed.gz', '.evidence.tsv')
                else:
                    outfile = infile.replace('.giggle.swap.intersect.bed', '.evidence.tsv')
                fut = ex.submit(intersect2evidence, infile, outfile, right_gene_col, sample_column, bgzip, burden)
                futures[fut] = (idx, outfile, col_evidence)
            for fut in as_completed(futures):
                idx, outfile, col = futures[fut]
                fut.result()
                df_evidence.at[idx, col] = outfile
    return df_evidence


def agg_evidence_by_category(outdir_g2f, outdir_agg, df_giggle_shards, outfile_suffix='-fusion_evidence.tsv', evidence_pattern='*.evidence.tsv'):
    '''
    aggregate giggle+bedtools evidence by category
    '''
    # validation
    assert os.path.exists(outdir_g2f), f"Input directory {outdir_g2f} does not exist"
    os.makedirs(outdir_agg, exist_ok=True)
    # parse categories
    categories = df_giggle_shards['category'].unique()
    k = len(categories)
    # verify no prior aggregation files exist
    for cat in categories:
        outfile_agg = os.path.join(outdir_agg, f"{cat}{outfile_suffix}")
        assert not os.path.exists(outfile_agg), f"Aggregation file {outfile_agg} already exists, please remove before running aggregation"
    # sequentially aggregate evidence within each category
    for i, cat in enumerate(categories):
        print(f"# category: {cat} ({i+1}/{k})")
        outdir_cat = os.path.join(outdir_g2f, cat)
        outfile_agg = os.path.join(outdir_agg, f"{cat}{outfile_suffix}")
        with open(outfile_agg, 'w') as f_out:
            # write header
            f_out.write(f'gene_left\tgene_right\treads_{cat}\tsamples_{cat}\n')
            # write all fusion evidence in category to single file
            evidence_files = glob.glob(os.path.join(outdir_cat, evidence_pattern))
            n = len(evidence_files)
            for j, path_evidence in enumerate(evidence_files):
                print(f"# processing file {j+1}/{n}: {path_evidence}")
                with open(path_evidence, 'r') as f_in:
                    # skip header
                    # gracefully handle empty evidence files
                    try:
                        next(f_in)
                    except StopIteration:
                        continue
                    # don't expect any comment lines, but just incase
                    for line in f_in:
                        if line.startswith('#'):
                            continue
                        else:
                            f_out.write(line)
    return outdir_agg
