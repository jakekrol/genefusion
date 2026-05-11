import pandas as pd
import polars as pl
import numpy as np
import tempfile
import subprocess
import shutil
import os
import shlex
from concurrent.futures import ThreadPoolExecutor, as_completed
from polymerization.io import *
# order

# 1. load and validate fusion set and bed
# 2. sort gene pairs in fusion set by their position in bed file
# 3. merge sorted fusion set with bed file data

# verify all fusion set genes are in bed file
def verify_fusion_set_in_bed(df_fusion, df_bed):
    '''
    raise assertion error if any gene in fusion set is not found in bed file
    otherwise, return None
    df_fusion: a 2 col tsv of gene names
    df_bed: a pandas dataframe of bed file, with at least 4 columns of chrom, start, end, and gene_name
    '''
    if isinstance(df_fusion, pl.DataFrame):
        fusion_cols = df_fusion.columns[:2]
        bed_col = 'gene_name'
        set_fusion_genes = set(df_fusion.select(fusion_cols).to_series(0).to_list()).union(
            set(df_fusion.select(fusion_cols).to_series(1).to_list())
        )
        set_bed_genes = set(df_bed.get_column(bed_col).to_list())
    else:
        set_fusion_genes = set(df_fusion.iloc[:, 0]).union(set(df_fusion.iloc[:, 1]))
        set_bed_genes = set(df_bed['gene_name'])
    test_subset = set_fusion_genes.issubset(set_bed_genes)
    if test_subset:
        return None
    else:
        missing_genes = set_fusion_genes - set_bed_genes
        raise AssertionError(f"The following genes in fusion set are not found in bed file: {missing_genes}")

def left_sort_fusion_set(df_fusion, df_bed):
    '''
    sort gene pairs by their position in bed file
    "smaller" gene is left, "larger" gene is right
    determined by start coordinate
    '''
    for i,row in df_fusion.iterrows():
        gene_x = row['gene_x']
        gene_y = row['gene_y']
        # read_bed already sorts the bed file by chromosome and start position
        # therefore, we can compare index in bed_file
        try:
            idx_x = df_bed.index[df_bed['gene_name'] == gene_x][0]
            idx_y = df_bed.index[df_bed['gene_name'] == gene_y][0]
        except IndexError:
            continue
        if idx_x < idx_y:
            continue
        else:
            # swap gene_x and gene_y
            df_fusion.loc[row.name, ['gene_x', 'gene_y']] = [gene_y, gene_x]
    cols = df_fusion.columns.tolist()
    df_fusion.columns = ['gene_left', 'gene_right'] + cols[2:]
    return df_fusion

def merge_fusion_set_with_bed(df_fusion, df_bed, merger='pandas'):
    '''
    merge fusion set with bed file data
    '''
    # check that all fusion set genes are in bed file
    verify_fusion_set_in_bed(df_fusion, df_bed)
    if merger == 'pandas':
        df_bed_to_merge = df_bed.copy()
        cols = [str(col) for col in df_bed_to_merge.columns.tolist()]
        left_cols = [col + '_left' for col in cols]
        right_cols = [col + '_right' for col in cols]
        df_bed_to_merge.columns = left_cols
        df_merged = pd.merge(df_fusion, df_bed_to_merge, left_on='gene_left', right_on='gene_name_left', how='left')
        df_bed_to_merge.columns = right_cols
        df_merged = pd.merge(df_merged, df_bed_to_merge, left_on='gene_right', right_on='gene_name_right', how='left', suffixes=('_left', '_right'))
        df_merged = df_merged.drop(columns=['gene_name_left', 'gene_name_right'])
    elif merger == 'polars':
        df_bed_left = df_bed.rename({col: f"{col}_left" for col in df_bed.columns})
        df_bed_right = df_bed.rename({col: f"{col}_right" for col in df_bed.columns})
        df_merged = df_fusion.join(df_bed_left, left_on='gene_left', right_on='gene_name_left', how='left')
        df_merged = df_merged.join(df_bed_right, left_on='gene_right', right_on='gene_name_right', how='left')
    return df_merged

def run_stix(argstring,left_gene, right_gene, outfile, timeout=60 * 60 * 2):
    # default timeout is 2 hours
    '''
    run stix with given arguments and write output to file
    '''
    stix = shutil.which('stix')
    if not stix:
        raise FileNotFoundError("stix command not found in PATH")

    cmd = [stix] + shlex.split(argstring)
    try:
        with open(outfile, 'w') as f:
            # write the genes in the header
            f.write(f"#gene_left={left_gene}\n")
            f.write(f"#gene_right={right_gene}\n")
            f.flush() # necessary to ensure header is written before suprocess runs
            print(f"Running command: {' '.join(cmd)}")
            subprocess.run(cmd, text=True, timeout=timeout, check=True, stdout=f, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running stix: {e}")
        print(f"Stderr: {e.stderr}")
        raise e
    except subprocess.TimeoutExpired as e:
        print(f"Stix command timed out after {timeout} seconds")
        raise e
    return True

def stix_shard_tbl2tmp(df_shard):
    '''
    write stix shard table to temporary tsv files for stix
    '''
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        # only write the index and peddb columns for stix
        df_shard_to_write = df_shard[['giggle_index', 'ped_db']]
        df_shard_to_write.to_csv(f.name, sep='\t', index=False, header=False)
        return f.name



def merge_fusion_set_bed2stix(
    df_merged,
    df_shard,
    outdir,
    outfile_prefix='',
    outfile_suffix='.stix',
    max_workers=4,
    gene_delim='--'
):
    '''
    run stix for every fusion in merged fusionset and bed file table
    '''
    # fusions are processed parallel
    # categories are processed sequential
    # same-category shards are processed in batches using stix -B shardfile

    records = df_merged.to_dict('records')
    max_workers=min(max_workers, len(records), os.cpu_count())

    for cat,group in df_shard.groupby('category'):
        # each category get its own directory
        outdir_cat = os.path.join(outdir, cat)
        os.makedirs(outdir_cat, exist_ok=True)
        alt_file_col = group['alt_file_col'].iloc[0]
        try:
            # the shardfile is used to aggregate queries across shards (still within same category)
            shardfile = stix_shard_tbl2tmp(group)
            # parallel stix queries
            with ThreadPoolExecutor(max_workers=max_workers) as ex:
                # for each fusion
                for r in records:
                    # query data
                    gene_left = r['gene_left']
                    gene_right = r['gene_right']

                    chromosome_left = r['chromosome_left']
                    chromosome_right = r['chromosome_right']
                    start_left = r['start_left']
                    start_right = r['start_right']
                    end_left = r['end_left']
                    end_right = r['end_right']
                    region_left = f"{chromosome_left}:{start_left}-{end_left}"
                    region_right = f"{chromosome_right}:{start_right}-{end_right}"

                    # construct stix command arguments
                    argstring = f"-l {region_left} -r {region_right} -B {shardfile} " \
                                f"-t FUSION -s 500 -c {alt_file_col}"
                    outfile = os.path.join(outdir_cat, f"{outfile_prefix}_{gene_left}{gene_delim}{gene_right}_{outfile_suffix}")
                    # run stix
                    ex.submit(run_stix, argstring, gene_left, gene_right, outfile)
                    
        finally:
            os.remove(shardfile)

            

def stix_fusion2evidence(df_stix_output):
    '''
    convert stix output for a fusion query into evidence format
    '''
    df_stix_output['total_evidence'] = df_stix_output['Pairend'] + df_stix_output['Split']
    reads = df_stix_output['total_evidence'].sum()
    samples = df_stix_output['total_evidence'].gt(0).sum()
    assert reads >= 0, "Read count should be non-negative"
    assert samples >= 0, "Sample count should be non-negative"
    return reads, samples


def agg_stix_evidence_by_category(outdir_s2f, outdir_agg, df_stix_shards, outfile_suffix="-fusion_evidence.tsv"):
    '''
    aggregate stix fusion evidence by category
    '''
    # validation
    assert os.path.exists(outdir_s2f), f"stix output directory does not exist: {outdir_s2f}"
    categories = df_stix_shards['category'].unique()
    k = len(categories)
    # gather all outfile_aggs to verify they do not exist prior to aggregation
    for cat in categories:
        outfile_agg = os.path.join(outdir_agg, f"{cat}{outfile_suffix}")
        assert not os.path.exists(outfile_agg), f"An aggregation file already exists: {outfile_agg}. Stopping to avoid double writes"
    # sequentially aggregate evidence for each category
    for i, cat in enumerate(categories):
        # location of stix files
        print("category:", cat)
        outdir_cat = os.path.join(outdir_s2f, cat)
        assert os.path.exists(outdir_cat), f"stix output category directory does not exist: {outdir_cat}"
        # file to write all fusion evidence to
        outfile_agg = os.path.join(outdir_agg, f"{cat}{outfile_suffix}")
        os.makedirs(outdir_agg, exist_ok=True)
        print(f"# aggregating evidence for category {i+1}/{len(categories)} named {cat} into file {outfile_agg}")
        with open(outfile_agg, 'w') as f_out:
            f_out.write(f"gene_left\tgene_right\treads_{cat}\tsamples_{cat}\n")
            # count reads and samples for each fusion
            n = len(os.listdir(outdir_cat))
            for j,filename in enumerate(os.listdir(outdir_cat)):
                filepath = os.path.join(outdir_cat, filename)
                print(f"# processing file {j+1}/{n} named {filename}")
                gene_left, gene_right, df_stix_output = read_stix_fusion_output(filepath)
                reads, samples = stix_fusion2evidence(df_stix_output)
                f_out.write(f"{gene_left}\t{gene_right}\t{reads}\t{samples}\n")

def fusion_tbl2_score_input(
    df_fusion,
    tumor_colmap,
    normal_colmap
):
    '''
    convert fusion table to score input format
    colmaps are dictionaries housing

    sub dictionaries for each subpopulation (tumor/normal)
    the sub dictionary requires the following keys: reads_col, samples_col, total_samples, and upper_bound
    '''
    M = df_fusion.shape[0]
    n_tumor_subpops = len(tumor_colmap.keys())
    n_normal_subpops = len(normal_colmap.keys())
    tumor_matrices = np.empty((M, n_tumor_subpops, 4), dtype=np.float64)
    normal_matrices = np.empty((M, n_normal_subpops, 4), dtype=np.float64)
    # col_tumor_blood_rna = df_fusion.iloc[:,3]
    for i,(key,subpop_dict) in enumerate(tumor_colmap.items()):
            reads_col = subpop_dict['reads_col']
            samples_col = subpop_dict['samples_col']
            total_samples = subpop_dict['total_samples']
            upper_bound = subpop_dict['upper_bound']
            tumor_matrices[:, i, 0] = df_fusion[reads_col]
            tumor_matrices[:, i, 1] = df_fusion[samples_col]
            tumor_matrices[:, i, 2] = total_samples
            tumor_matrices[:, i, 3] = upper_bound
    for i,(key,subpop_dict) in enumerate(normal_colmap.items()):
            reads_col = subpop_dict['reads_col']
            samples_col = subpop_dict['samples_col']
            total_samples = subpop_dict['total_samples']
            upper_bound = subpop_dict['upper_bound']
            normal_matrices[:, i, 0] = df_fusion[reads_col]
            normal_matrices[:, i, 1] = df_fusion[samples_col]
            normal_matrices[:, i, 2] = total_samples
            normal_matrices[:, i, 3] = upper_bound
    return tumor_matrices, normal_matrices


    

