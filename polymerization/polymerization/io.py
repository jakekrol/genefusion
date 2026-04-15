import pandas as pd
import os

def validate_fusion_set(df_fusion):
    '''
    raise assertion error if duplicates or selfies are found in input fusion set
    otherwise, return None
    df_fusion: a 2 col tsv of gene names
    '''
    # check for duplicates
    # keep=False will also mark the first occurrence of duplicates
    series_bool = df_fusion.duplicated(keep=False)
    assert not(any(series_bool))
    # check for selfies
    series_bool = df_fusion.iloc[:, 0] == df_fusion.iloc[:, 1]
    assert not(any(series_bool))
    # check column names
    assert len(df_fusion.columns) >= 2
    required_columns = ['gene_x', 'gene_y']
    for col in required_columns:
        assert col in df_fusion.columns
    return None

def read_fusion_set(path, header=None, sep='\t'):
    '''
    read in fusion set from tsv file and return as pandas dataframe
    path_fusion_set: path to tsv file of fusion set, with 2 columns of gene names
    '''
    df_fusion = pd.read_csv(path, sep=sep, header=header)
    df_fusion.columns=['gene_x', 'gene_y']
    return df_fusion

def validate_bed(df_bed):
    '''
    raise assertion error if bed file is not in correct format
    otherwise, return None
    df_bed: a pandas dataframe of bed file, with at least 4 columns of chromosome, start, end, and gene_name
    '''
    # check for duplicate rows
    series_bool = df_bed.duplicated(keep=False)
    assert not(any(series_bool))
    # check for essential columns
    assert len(df_bed.columns) >= 4
    required_cols = ['chromosome', 'start', 'end', 'gene_name']
    for col in required_cols:
        assert col in df_bed.columns
    # check datatype of start and end columns
    assert pd.api.types.is_integer_dtype(df_bed['start'])
    assert pd.api.types.is_integer_dtype(df_bed['end'])
    # start should be less than end
    series_bool = df_bed['start'] <= df_bed['end']
    assert all(series_bool)
    return None

def read_bed(path, gene_col_idx, header=None, sep='\t'):
    '''
    read in bed file and return as pandas dataframe
    path_bed: path to bed file, with 4 columns of chromosome, start, end, and gene_name
    gene_col_idx: index (0-start) of the column containing gene names
    '''
    df_bed = pd.read_csv(path, sep=sep, header=header)
    cols = df_bed.columns.tolist()
    cols = ['chromosome', 'start', 'end'] + cols[3:]
    cols[gene_col_idx] = 'gene_name'
    df_bed.columns = cols
    # sort by chromosome and start position
    df_bed['chromosome'] = df_bed['chromosome'].astype(str)
    df_bed = df_bed.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    validate_bed(df_bed)
    return df_bed

def validate_stix_shardfile(df_stix_shards):
    '''
    raise assertion error if stix shard file is not in correct format
    otherwise, return None
    df_stix_shards: a pandas dataframe of stix shard file, with 3 columns of giggle_index, ped_db, and category
    '''
    # check for duplicate rows
    series_bool = df_stix_shards.duplicated(keep=False)
    assert not(any(series_bool))
    # check for essential columns
    assert len(df_stix_shards.columns) >= 3
    required_cols = ['giggle_index', 'ped_db', 'category']
    for col in required_cols:
        assert col in df_stix_shards.columns
    # verify all paths exist
    for i in df_stix_shards.index:
        try:
            path_giggle_index = df_stix_shards.loc[i, 'giggle_index']
            assert os.path.exists(path_giggle_index)
        except AssertionError:
            raise AssertionError(f"Path to giggle index not exist for row {i} of shard file")
        try:
            path_ped_db = df_stix_shards.loc[i, 'ped_db']
            assert os.path.exists(path_ped_db)
        except AssertionError:
            raise AssertionError(f"Path to ped db not exist for row {i} of shard file")
    return None

def read_stix_shardfile(path, sep='\t', header=None):
    '''
    read in stix shard file and return as pandas dataframe
    alt_file_col is 1-indexed check the ped file for this
    '''
    df_stix_shards = pd.read_csv(path, sep=sep, header=header)
    if not header:
        # assign column names if not provided
        df_stix_shards.columns = ['giggle_index', 'ped_db', 'alt_file_col', 'category']
    validate_stix_shardfile(df_stix_shards)
    return df_stix_shards

def validate_stix_fusion_output(df_stix_output):
    '''
    check for necessary columns in stix output file
    '''
    required_cols = ['Giggle_File_Id', 'Pairend', 'Split']
    for col in required_cols:
        assert col in df_stix_output.columns, f"Missing required column in stix output file: {col}"
    # check column types
    # assert pd.api.types.is_integer_dtype(df_stix_output['Pairend']), "Pairend column should be of integer type"
    # assert pd.api.types.is_integer_dtype(df_stix_output['Split']), "Split column should be of integer type"

def read_stix_fusion_output(path,misc_header_line=2, sep='\t'):
    '''
    read in stix output file and return as pandas dataframe
    misc_header line: 0-indexed line number of the stix default header
    sep: delimiter
    returns a tuple of (gene_left, gene_right, df_stix_output)
    '''
    # get the query genes from header
    gene_left,gene_right = None, None
    with open(path, 'r') as f:
        for line in f:
            # stop reading once we exhaust header lines
            if not line.startswith('#'):
                break
            if line.startswith('#gene_left='):
                gene_left = line.strip().split('=')[1]
            if line.startswith('#gene_right='):
                gene_right = line.strip().split('=')[1]

    df_stix_output = pd.read_csv(path, sep=sep, skiprows=[misc_header_line], comment='#')
    validate_stix_fusion_output(df_stix_output)
    return (gene_left, gene_right, df_stix_output)

