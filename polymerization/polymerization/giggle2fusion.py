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
import time

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


def _is_bad_chrom_token(token):
    token = token.lower()
    return (
        token == "0"
        or token.startswith("hs")
        or token.startswith("gl")
        or token.startswith("nc")
        or token.startswith("mt")
        or token.startswith("-1")
        or token.startswith("*")
    )

def merge_fusion_set_bed2giggle(
    df_merged,
    df_shard,
    outdir,
    outfile_prefix='',
    max_workers=4,
    gene_delim='--',
    timeout=60 * 60 * 2,
    bgzip=True,
    outfile_column_giggle_prefix = 'outfile_giggle'
):
    # similar to stix2fusion.merge_fusion_set_bed2stix() but for giggle instead of stix

    # only left-genes are searched
    # significant reduction in complexity
    right_gene_cols = [df_merged.columns.get_loc(col) for col in df_merged.columns if 'right' in col]
    df_merged.drop(columns=[df_merged.columns[i] for i in right_gene_cols], inplace=True)
    df_merged = df_merged.drop_duplicates(subset='gene_left', keep='first')
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
                        outdir_cat, f"{outfile_prefix}{gene_left}.giggle"
                    )
                    if bgzip and not outfile.endswith('.gz'):
                        outfile += '.gz'
                    # run giggle
                    futures.append(ex.submit(run_giggle, argstring, gene_left, outfile, timeout, bgzip, append=True))
                for future in as_completed(futures):
                    try:
                        res = future.result()
                    except Exception as e:
                        print(f"Error in giggle search: {e}")
                        continue
                    if res is not None:
                        df_merged.loc[df_merged['gene_left'] == gene, col_name] = outfile
                    else:
                        df_merged.loc[df_merged['gene_left'] == gene, col_name] = pd.NA
    # includes giggle outfile column
    df_giggle = df_merged.copy()
    return df_giggle

def clean_excord(path_excord, outfile, bgzip=False):
    # check if file is empty or dne
    if not os.path.exists(path_excord) or os.path.getsize(path_excord) == 0:
        print(f"Warning: Excord file {path_excord} does not exist or is empty, skipping cleaning")
        return None
    input_open = gzip.open if bgzip else open
    output_open = pysam.BGZFile if bgzip else open
    input_mode = 'rt' if bgzip else 'r'
    output_mode = 'wb' if bgzip else 'w'

    with input_open(path_excord, input_mode) as f_in:
        with output_open(outfile, output_mode) as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line.encode() if bgzip else line)
                    continue

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 10:
                    continue

                left_chrom = fields[0]
                left_start = fields[1]
                right_chrom = fields[4]
                right_start = fields[5]

                if _is_bad_chrom_token(left_chrom) or _is_bad_chrom_token(right_chrom):
                    continue
                if (left_chrom.lower() == "0") and (left_start == "0"):
                    continue
                if (right_chrom.lower() == "0") and (right_start == "0"):
                    continue

                f_out.write(line.encode() if bgzip else line)
    return outfile

def giggle2clean(
    df_giggle,
    df_shard,
    bgzip=False,
    outfile_column_giggle_prefix='outfile_giggle',
    outfile_column_clean_prefix='outfile_clean',
    max_workers=4
):
    df_clean = df_giggle.copy()
    for cat in df_shard['category'].unique():
        col_giggle = f"{outfile_column_giggle_prefix}_{cat}"
        col_clean = f"{outfile_column_clean_prefix}_{cat}"
        df_clean[col_clean] = '' # initialize clean outfile column for this category
        # parallel clean excord output
        workers = min(max_workers, os.cpu_count() or 1, len(df_giggle))
        futures = {}
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for idx, row in df_giggle.iterrows():
                infile = row[f"{col_giggle}"]
                if bgzip:
                    outfile = infile.replace('.giggle.gz', '.giggle.clean.gz')
                else:
                    outfile = infile.replace('.giggle', '.giggle.clean')
                fut = ex.submit(clean_excord, infile, outfile, bgzip)
                futures[fut] = (idx, outfile, col_clean)
            for fut in as_completed(futures):
                idx, outfile, col = futures[fut]
                try:
                    res = fut.result()
                # worker fails
                except Exception as e:
                    print(f"giggle2clean: worker failed for idx={idx}, infile={df_giggle.at[idx, col_giggle]}: {e}")
                    df_clean.at[idx, col] = pd.NA #
                    continue
                # worker succeeds, but clean_excord returns None (e.g. empty or dne file)
                if res is None:
                    df_clean.at[idx, col] = pd.NA
                else: 
                    df_clean.at[idx, col] = outfile
    return df_clean

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
    if not os.path.exists(path_gigglefile) or os.path.getsize(path_gigglefile) == 0:
        return None
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
    return outfile

def clean2swap(
    df_clean,
    df_shard,
    bgzip=False,
    outfile_column_clean_prefix='outfile_clean',
    outfile_column_swap_prefix='outfile_swap',
    max_workers=4
):
    df_swap = df_clean.copy()
    for cat in df_shard['category'].unique():
        col_clean = f"{outfile_column_clean_prefix}_{cat}"
        col_swap = f"{outfile_column_swap_prefix}_{cat}"
        df_swap[col_swap] = '' # initialize swap outfile column for this category
        # parallel swap left and right intervals in giggle output
        workers = min(max_workers, os.cpu_count() or 1, len(df_clean))
        futures = {}
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for idx, row in df_clean.iterrows():
                infile = row[col_clean]
                outfile = infile.replace('.giggle.clean', '.giggle.clean.swap')
                fut = ex.submit(swap_intervals, infile, outfile, bgzip)
                futures[fut] = (idx, outfile, col_swap)
            for fut in as_completed(futures):
                idx, outfile, col = futures[fut]
                try:
                    res = fut.result()
                except Exception as e:
                    print(f"clean2swap: worker failed for idx={idx}, infile={df_clean.at[idx, col_clean]}: {e}")
                    df_swap.at[idx, col] = pd.NA
                    continue
                if res is None:
                    df_swap.at[idx, col] = pd.NA
                else:
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
    if not os.path.exists(path_input) or os.path.getsize(path_input) == 0:
        return None
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
    return outfile

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
                outfile = infile.replace('.giggle.clean.swap', '.giggle.clean.swap.intersect.bed')
                fut = ex.submit(bedtools_intersect, infile, path_bedfile, outfile, gene_col_idx, bgzip, bedtools_bin)
                futures[fut] = (idx, outfile, col_intersect)
            for fut in as_completed(futures):
                idx, outfile, col = futures[fut]
                try:
                    res = fut.result()
                except Exception as e:
                    print(f"swap2intersect: worker failed for idx={idx}, infile={df_swap.at[idx, col_swap]}: {e}")
                    df_intersect.at[idx, col] = pd.NA
                    continue
                if res is None:
                    df_intersect.at[idx, col] = pd.NA
                else:
                    df_intersect.at[idx, col] = outfile
    return df_intersect

    
def intersect2evidence(
    path_intersect,
    outfile,
    right_gene_col=3,
    sample_column=14,
    bgzip=False,
    burden=False,
    sample_clean_func=None,
):
    if not os.path.exists(path_intersect) or os.path.getsize(path_intersect) == 0:
        return None
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
    # set default sample id cleaning function
    if sample_clean_func is None:
        sample_clean_func = lambda x: os.path.basename(x)

    # handles bgzipped or not automatically
    df = pd.read_csv(path_intersect, sep='\t', header=None, comment='#', usecols=[right_gene_col, sample_column])
    df.columns=['gene_right', 'sample_id']
    # read counts
    right_gene_counts = df.groupby('gene_right').size().reset_index(name='reads')
    # sample counts
    # apply sample id cleaning function if provided
    df['sample_id'] = df['sample_id'].apply(sample_clean_func)
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
            f.write("#gene_left={}\n".format(gene_left))
            f.write(f"{total_burden}\n")
    # write output
    df_evidence.to_csv(outfile, sep='\t', index=False, compression=None)
    return outfile

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
                infile = str(row[col_intersect])
                outfile = infile.replace('.giggle.clean.swap.intersect.bed.gz', '.evidence.tsv')
                fut = ex.submit(intersect2evidence, infile, outfile, right_gene_col, sample_column, bgzip, burden)
                futures[fut] = (idx, outfile, col_evidence)
            for fut in as_completed(futures):
                idx, outfile, col = futures[fut]
                try:
                    res = fut.result()
                except Exception as e:
                    print(f"intersect2evidence: worker failed for idx={idx}, infile={df_intersect.at[idx, col_intersect]}: {e}")
                    df_evidence.at[idx, col] = pd.NA
                    continue
                if res is None:
                    df_evidence.at[idx, col] = pd.NA
                else:
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

def giggle2fusion(
    df_merged,
    df_shard,
    outdir,
    logdir,
    path_bedfile,
    gene_col_idx=3,
    sample_clean_func= lambda x: os.path.basename(x), # default takes basename of sample id
    evidence_right_gene_col=3,
    evidence_sample_col=14,
    outfile_prefix='',
    max_workers=4,
    timeout=60 * 60 * 2,
    bgzip=True,
    bedtools_bin=None,
    burden=True,
    verbose=False
):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)
    with open(os.path.join(logdir, 'giggle2fusion_parameters.txt'), 'w') as f:
        f.write(f"outdir: {outdir}\n")
        f.write(f"logdir: {logdir}\n")
        f.write(f"path_bedfile: {path_bedfile}\n")
        f.write(f"gene_col_idx: {gene_col_idx}\n")
        f.write(f"sample_clean_func: {sample_clean_func}\n")
        f.write(f"evidence_right_gene_col: {evidence_right_gene_col}\n")
        f.write(f"evidence_sample_col: {evidence_sample_col}\n")
        f.write(f"outfile_prefix: {outfile_prefix}\n")
        f.write(f"max_workers: {max_workers}\n")
        f.write(f"timeout: {timeout}\n")
        f.write(f"bgzip: {bgzip}\n")
        f.write(f"bedtools_bin: {bedtools_bin}\n")
        f.write(f"burden: {burden}\n")
        f.write(f"verbose: {verbose}\n")
    with open(os.path.join(logdir, 'NA_counter.tsv'), 'w') as f:
        f.write(f"step\tNA_count\n")

    # giggle search
    if verbose:
        print("# running giggle search")
    t_0= time.time()
    df_giggle = merge_fusion_set_bed2giggle(
        df_merged, df_shard, outdir, outfile_prefix, max_workers, gene_delim='--', timeout=timeout, bgzip=bgzip
    )
    with open (os.path.join(logdir, 'giggle_time.txt'), 'w') as f:
        f.write(f"{time.time() - t_0}\n")
    df_giggle.to_csv(os.path.join(logdir, 'giggle_output.tsv'), sep='\t', index=False)
    with open(os.path.join(logdir, 'NA_counter.tsv'), 'a') as f:
        f.write(f"giggle_output\t{df_giggle.isnull().sum().sum()}\n")

    # clean
    if verbose:
        print("# cleaning giggle output")
    t_0=    time.time()
    df_clean = giggle2clean(
        df_giggle, df_shard, bgzip=bgzip, outfile_column_giggle_prefix='outfile_giggle', outfile_column_clean_prefix='outfile_clean', max_workers=max_workers
    )
    with open(os.path.join(logdir, 'giggle_clean_time.txt'), 'w') as f:
        f.write(f"{time.time() - t_0}\n")
    df_clean.to_csv(os.path.join(logdir, 'giggle_clean_output.tsv'), sep='\t', index=False)
    with open(os.path.join(logdir, 'NA_counter.tsv'), 'a') as f:
        f.write(f"giggle_clean_output\t{df_clean.isnull().sum().sum()}\n")

    # swap
    if verbose:
        print("# swapping intervals in giggle output")
    t_0= time.time()
    df_swap = clean2swap(
        df_clean,
        df_shard,
        bgzip=bgzip,
        outfile_column_clean_prefix='outfile_clean',
        outfile_column_swap_prefix='outfile_swap',
        max_workers=max_workers
    )
    with open(os.path.join(logdir, 'giggle_swap_time.txt'), 'w') as f:
        f.write(f"{time.time() - t_0}\n")
    df_swap.to_csv(os.path.join(logdir, 'giggle_swap_output.tsv'), sep='\t', index=False)
    with open(os.path.join(logdir, 'NA_counter.tsv'), 'a') as f:
        f.write(f"giggle_swap_output\t{df_swap.isnull().sum().sum()}\n")

    # bedtools intersect
    if verbose:
        print("# intersecting swapped intervals with bed file")
    t_0= time.time()
    df_intersect = swap2intersect(
        df_swap, df_shard, path_bedfile, gene_col_idx, outfile_column_swap_prefix='outfile_swap', outfile_column_intersect_prefix='outfile_intersect', bgzip=bgzip, bedtools_bin=bedtools_bin, max_workers=max_workers
    )
    with open(os.path.join(logdir, 'giggle_intersect_time.txt'), 'w') as f:
        f.write(f"{time.time() - t_0}\n")
    df_intersect.to_csv(os.path.join(logdir, 'giggle_intersect_output.tsv'), sep='\t', index=False)
    with open(os.path.join(logdir, 'NA_counter.tsv'), 'a') as f:
        f.write(f"giggle_intersect_output\t{df_intersect.isnull().sum().sum()}\n")

    # intersect to evidence
    if verbose:
        print("# converting intersect files to evidence files")
    t_0= time.time()
    df_evidence = df_intersect2df_evidence(
        df_intersect, df_shard, evidence_right_gene_col, evidence_sample_col,
        outfile_column_intersect_prefix='outfile_intersect', outfile_column_evidence_prefix='outfile_evidence',
        bgzip=bgzip, burden=burden, max_workers=max_workers
    )
    with open(os.path.join(logdir, 'giggle_evidence_time.txt'), 'w') as f:
        f.write(f"{time.time() - t_0}\n")
    df_evidence.to_csv(os.path.join(logdir, 'giggle_evidence_output.tsv'), sep='\t', index=False)
    with open(os.path.join(logdir, 'NA_counter.tsv'), 'a') as f:
        f.write(f"giggle_evidence_output\t{df_evidence.isnull().sum().sum()}\n")

    return df_evidence
    