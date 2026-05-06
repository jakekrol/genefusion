# see stix2fusion for preprocessing steps
import shlex
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from polymerization.io import *
import os
import pandas as pd
import numpy as np

def run_giggle(argstring, left_gene, right_gene, outfile, timeout = 60 * 60 * 2):
    '''
    run giggle with given arguments and return the output as a pandas dataframe
    argstring: string of arguments to pass to giggle, should include {left_gene} and {right_gene} as placeholders for the gene names
    left_gene: name of left gene
    right_gene: name of right gene
    outfile: path to output file to write giggle results to
    timeout: time limit for giggle to run in seconds, default is 2 hours
    returns a pandas dataframe of giggle results
    '''
    giggle = shutil.which('giggle')
    if not giggle:
        raise FileNotFoundError("giggle command not found in PATH")
    cmd = [giggle] + shlex.split(argstring)
    # run the command using subprocess.run with a timeout
    try:
        with open(outfile, 'w') as f:
            # write genes in the header
            f.write(f"#gene_left={left_gene}\n")
            f.write(f"#gene_right={right_gene}\n")
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
    return True

def merge_fusion_set_bed2giggle(
    df_merged,
    df_shard,
    outdir,
    outfile_prefix='',
    outfile_suffix='.giggle',
    max_workers=4,
    gene_delim='--'
):
    # similar to stix2fusion.merge_fusion_set_bed2stix() but for giggle instead of stix
    records = df_merged.to_dict('records') # list of row dictionaries
    max_workers=min(max_workers, len(records), os.cpu_count())

    for cat,group in df_shard.groupby('category'):
        # each category get its own directory
        outdir_cat = os.path.join(outdir, cat)
        os.makedirs(outdir_cat, exist_ok=True)
        giggle_index = group['giggle_index'].iloc[0]
        # parallel giggle queries
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            # for each fusion
            for r in records:
                # query data
                gene_left = r['gene_left']
                gene_right = r['gene_right']

                chromosome_left = r['chromosome_left']
                # chromosome_right = r['chromosome_right']
                start_left = r['start_left']
                # start_right = r['start_right']
                end_left = r['end_left']
                # end_right = r['end_right']
                region_left = f"{chromosome_left}:{start_left}-{end_left}"
                # region_right = f"{chromosome_right}:{start_right}-{end_right}"

                # construct giggle command arguments
                argstring = f"search -i {giggle_index} -r {region_left} -v"
                outfile = os.path.join(
                    outdir_cat, f"{outfile_prefix}_{gene_left}{gene_delim}{gene_right}_{outfile_suffix}"

                )
                # run giggle
                ex.submit(run_giggle, argstring, gene_left, gene_right, outfile)