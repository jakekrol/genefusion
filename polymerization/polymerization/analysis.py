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