import pandas as pd
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
    set_fusion_genes = set(df_fusion.iloc[:, 0]).union(set(df_fusion.iloc[:, 1]))
    set_bed_genes = set(df_bed['gene_name'])
    assert set_fusion_genes.issubset(set_bed_genes)
    return None

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
        idx_x = df_bed.index[df_bed['gene_name'] == gene_x][0]
        idx_y = df_bed.index[df_bed['gene_name'] == gene_y][0]
        if idx_x < idx_y:
            continue
        else:
            # swap gene_x and gene_y
            df_fusion.loc[row.name, ['gene_x', 'gene_y']] = [gene_y, gene_x]
    df_fusion.columns = ['gene_left', 'gene_right']
    return df_fusion

def merge_fusion_set_with_bed(df_fusion, df_bed):
    '''
    merge fusion set with bed file data
    '''
    df_merged = pd.merge(df_fusion, df_bed, left_on='gene_left', right_on='gene_name', how='left')
    df_merged = pd.merge(df_merged, df_bed, left_on='gene_right', right_on='gene_name', how='left', suffixes=('_left', '_right'))
    df_merged = df_merged.drop(columns=['gene_name_left', 'gene_name_right'])
    return df_merged
