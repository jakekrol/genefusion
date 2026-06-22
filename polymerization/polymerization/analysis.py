import pandas as pd

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
	T= np.zeros((bins_per_dim, bins_per_dim), dtype=np.int32)
	N= np.zeros((bins_per_dim, bins_per_dim), dtype=np.int32)
	for df_bp in list_df_bp_tumor:
		T += bin_breakpoints(df_bp, x_edges, y_edges)
	for df_bp in list_df_bp_normal:
		N += bin_breakpoints(df_bp, x_edges, y_edges)
	# breakpoint aware evidence is the minimum of T and N in each bin
	N_bp_aware = np.minimum(N, T)
	evidence = N_bp_aware.sum()
	return evidence