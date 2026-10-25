import gzip
import os
import pandas as pd
import numpy as np
import time
from collections import defaultdict


def intersect_file2sample_counts(path: str, gene_right_query: str, bgzip: bool = False):
    assert os.path.exists(path)
    ### relevant columns (0-indexed)
    COL_GENE_RIGHT=3
    COL_SAMPLE=14
    sample_counts=defaultdict(int)
    with gzip.open(path, "rt") as f:
        for i,line in enumerate(f.readlines()):
            if line.startswith("#"):
                continue
            columns = line.strip().split()
            gene_right = columns[COL_GENE_RIGHT]
            # ignore rows that do not match query
            if gene_right != gene_right_query:
                continue
            sample = columns[COL_SAMPLE]
            sample_counts[sample] += 1
    return sample_counts

def intersect_file2data(path: str, gene_right_query: str, group_by_sample: bool = False, bgzip=False,timeout=180):
    '''
    returns pd.DataFrame if successful; returns None if timeout or exception is caught
    '''
    assert os.path.exists(path)
    ### relevant columns (0-indexed)
    COL_GENE_RIGHT=3
    COL_SAMPLE=14
    # read info
    COL_CHROM_LEFT=9
    COL_START_LEFT=10
    COL_END_LEFT=11
    COL_STRAND_LEFT=12
    COL_CHROM_RIGHT=5
    COL_START_RIGHT=6
    COL_END_RIGHT=7
    COL_STRAND_RIGHT=8
    row_data = []
    try:
        opener = gzip.open if bgzip else open
        with opener(path, "rt") as f:
            columns=None
            t_0=time.time()
            for i,line in enumerate(f.readlines()):
                if (t_0 - time.time()) > timeout:
                    print(f"# intersect_file2data timed out for query: {path}")
                    return None
                if line.startswith("#"):
                    # gene_left=line.strip().split("=")[1]
                    continue
                columns = line.strip().split()
                gene_right = columns[COL_GENE_RIGHT]
                # ignore rows that do not match query
                if gene_right != gene_right_query:
                    continue
                sample = columns[COL_SAMPLE]
                chrom_left = columns[COL_CHROM_LEFT]
                start_left = int(columns[COL_START_LEFT])
                end_left = int(columns[COL_END_LEFT])
                strand_left = int(columns[COL_STRAND_LEFT])
                chrom_right = columns[COL_CHROM_RIGHT]
                start_right = int(columns[COL_START_RIGHT])
                end_right = int(columns[COL_END_RIGHT])
                strand_right = int(columns[COL_STRAND_RIGHT])
                row_data.append(
                    (sample, chrom_left, start_left, end_left, strand_left,
                     chrom_right, start_right, end_right, strand_right)
                )
        df = pd.DataFrame(
            row_data,
            columns = [
                'sample', 'chromosome_left', 'start_left', 'end_left', 'strand_left',
                'chromosome_right', 'start_right', 'end_right', 'strand_right'
            ]
        )
        if group_by_sample:
            df = df.groupby('sample').agg(
                start_left=('start_left', 'min'),
                end_right=('end_right', 'max')
            ).reset_index()
            # add back chromosomes, should be same in loop so using last iteration value is fine
            df['chromosome_left'] = chrom_left
            df['chromosome_right'] = chrom_right
        return df
    except Exception as e:
        print(f'# error in intersect_file2breakpoints, query={path}:', e)
        return None

# fast alternative to the read_g2f_intersect(); intersect2breakpoints() approach
def intersect_file2breakpoints(path: str, gene_right_query: str, group_by_sample: bool = False, bgzip=False,timeout=180):
    '''
    returns pd.DataFrame if successful; returns None if timeout or exception is caught
    '''
    assert os.path.exists(path)
    ### relevant columns (0-indexed)
    COL_GENE_RIGHT=3
    COL_SAMPLE=14
    # read info
    COL_CHROM_LEFT=9
    COL_START_LEFT=10
    COL_END_LEFT=11
    COL_CHROM_RIGHT=5
    COL_START_RIGHT=6
    COL_END_RIGHT=7
    row_data = []
    try:
        opener = gzip.open if bgzip else open
        with opener(path, "rt") as f:
            sample_breakpoints=defaultdict(dict)
            columns=None
            t_0=time.time()
            for i,line in enumerate(f.readlines()):
                if (t_0 - time.time()) > timeout:
                    print(f"# intersect_file2breakpoint timed out for query: {path}")
                    return None
                if line.startswith("#"):
                    continue
                columns = line.strip().split()
                gene_right = columns[COL_GENE_RIGHT]
                # ignore rows that do not match query
                if gene_right != gene_right_query:
                    continue
                sample = columns[COL_SAMPLE]
                chrom_left = columns[COL_CHROM_LEFT]
                start_left = int(columns[COL_START_LEFT])
                end_left = int(columns[COL_END_LEFT])
                chrom_right = columns[COL_CHROM_RIGHT]
                start_right = int(columns[COL_START_RIGHT])
                end_right = int(columns[COL_END_RIGHT])
                row_data.append(
                    (sample, chrom_left, start_left, end_left, chrom_right, start_right, end_right)
                )
        df = pd.DataFrame(
            row_data,
            columns = [
                'sample', 'chromosome_left', 'start_left', 'end_left',
                'chromosome_right', 'start_right', 'end_right'
            ]
        )
        if group_by_sample:
            df = df.groupby('sample').agg(
                start_left=('start_left', 'min'),
                end_right=('end_right', 'max')
            ).reset_index()
            # add back chromosomes, should be same in loop so using last iteration value is fine
            df['chromosome_left'] = chrom_left
            df['chromosome_right'] = chrom_right
        return df
    except Exception as e:
        print(f'# error in intersect_file2breakpoints, query={path}:', e)
        return None

def estimate_breakpoint(df_intersect: pd.DataFrame):
    # make sure to use group_by_sample=True for this
    assert "start_left" in df_intersect.columns
    assert "end_right" in df_intersect.columns
    left_chrom = df_intersect.loc[0,'chromosome_left']
    right_chrom = df_intersect.loc[0,'chromosome_right']
    bp_left = int(df_intersect['start_left'].mean())
    bp_right = int(df_intersect['end_right'].mean())
    bp_left = f"{left_chrom}:{bp_left}"
    bp_right = f'{right_chrom}:{bp_right}'
    return bp_left, bp_right
    
        

def intersect2chromosome_count(df_intersect):
    '''
    df_intersect: dataframe from read_g2f_intersect
    return count of intra and inter chromosomal fusion supporting reads
    '''
    left_gene_chrom = str(df_intersect['left_chromosome'].values[0])
    df_intersect['right_chromosome'] = df_intersect['right_chromosome'].astype(str)
    intra_count = (df_intersect['right_chromosome'] == left_gene_chrom).sum()
    inter_count = (df_intersect['right_chromosome'] != left_gene_chrom).sum()
    return {'intra': intra_count, 'inter': inter_count}


def intersect2breakpoints(df_intersect, gene_right, group_by_sample=False):
    '''
    extract left and right breakpoint positions from g2f intersect dataframe
    '''
    df_intersect = df_intersect[df_intersect['gene_right_name'] == gene_right]
    required_cols = [
        'left_start',
        'right_end'
    ]
    if group_by_sample:
        required_cols.append('sample')
    assert all(col in df_intersect.columns for col in required_cols), f"DataFrame must contain columns: {required_cols}"
    # optionally group by sample and take min left_start and max right_end for each sample
    if group_by_sample:
        df_breakpoints = df_intersect.groupby('sample').agg(
            left_breakpoint=('left_start', 'min'),
            right_breakpoint=('right_end', 'max')
        ).reset_index()
    else:
        df_breakpoints = df_intersect.copy()
        df_breakpoints['left_breakpoint'] = df_breakpoints['left_start']
        df_breakpoints['right_breakpoint'] = df_breakpoints['right_end']
    # move breakpoints to front
    cols = df_breakpoints.columns.tolist()
    cols.remove('left_breakpoint')
    cols.remove('right_breakpoint')
    cols = ['left_breakpoint', 'right_breakpoint'] + cols
    df_breakpoints = df_breakpoints[cols]
    return df_breakpoints

def bin_breakpoints(df_bp, x_edges, y_edges):
	'''
	bin breakpoints into 2D histogram
	'''
	hist, _, _ = np.histogram2d(df_bp['left_breakpoint'], df_bp['right_breakpoint'], bins=(x_edges, y_edges))
	# swap axes to have ERG on x and TMPRSS2 on y
	# np.histogram2d array puts x on rows and y on columns, so need to transpose to have ERG on x and TMPRSS2 on y
	hist = hist.T
	return hist

# for reads use read_g2f_intersect(..., group_by_sample=False)
# for samples use read_g2f_intersect(..., group_by_sample=True)
# as inputs
def count_breakpoint_aware_normal_evidence(list_df_bp_tumor,list_df_bp_normal, x_edges, y_edges, bins_per_dim=10):
	'''
	list_df_bp_tumor: list of dataframes of tumor breakpoints
	list_df_bp_normal: list of dataframes of normal breakpoints
	x_edges: edges for binning left breakpoints
	y_edges: edges for binning right breakpoints
	bins_per_dim: number of bins per dimension for histogram
	return number of breakpoint aware reads and samples in normal that are also in tumor
	'''
	T= np.zeros((bins_per_dim, bins_per_dim))
	N= np.zeros((bins_per_dim, bins_per_dim))
	for df_bp in list_df_bp_tumor:
		T += bin_breakpoints(df_bp, x_edges, y_edges)
	for df_bp in list_df_bp_normal:
		N += bin_breakpoints(df_bp, x_edges, y_edges)
	# breakpoint aware evidence is the minimum of T and N in each bin
	N_bp_aware = np.minimum(N, T)
	evidence = N_bp_aware.sum()
	return evidence