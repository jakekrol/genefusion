import pandas as pd
import pytest
from polymerization.io import *


# fusion set
def test_validate_fusion_set_fail():
    x = pd.DataFrame(
        {
            'gene_x': ['dup', 'dup', 'self'], 'gene_y': ['dup_partner', 'dup_partner', 'self']
        }
    )

    # only pass if error is raised
    with pytest.raises(AssertionError):
        validate_fusion_set(x)
    


def test_validate_fusion_set_pass():
    x = pd.DataFrame(
        {
            'gene_x': ['a', 'a', 'b'], 'gene_y': ['b', 'c', 'c']
        }
    )
    validate_fusion_set(x)


def test_validate_bed_fail():
    # missing gene_name column
    x = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr2'], 'start': [100, 150], 'end': [149, 300] 
        }
    )
    with pytest.raises(AssertionError):
        validate_bed(x)
    # duplicates
    x = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr1'], 'start': [100, 100], 'end': [149, 149], 'gene_name': ['geneA', 'geneA']
        }
    )
    with pytest.raises(AssertionError):
        validate_bed(x)
    # start > end
    x = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr2'], 'start': [150, 300], 'end': [149, 250], 'gene_name': ['geneA', 'geneB']
        }
    )
    with pytest.raises(AssertionError):
        validate_bed(x)

def test_validate_bed_pass():
    x = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr2'], 'start': [100, 150], 'end': [200, 300], 'gene_name': ['geneA', 'geneB']
        }
    )
    validate_bed(x)

def test_validate_stix_shardfile_pass(tmp_path):
    index1 = tmp_path / 'index1'
    index2 = tmp_path / 'index2'
    ped1 = tmp_path / 'db1'
    ped2 = tmp_path / 'db2'
    index1.mkdir()
    index2.mkdir()
    ped1.mkdir()
    ped2.mkdir()

    x = pd.DataFrame(
        {
            'giggle_index': [str(index1), str(index2)], 'ped_db': [str(ped1), str(ped2)], 'category': ['cat1', 'cat2']
        }
    )
    validate_stix_shardfile(x)

def test_validate_stix_shardfile_fail(tmp_path):
    index1 = tmp_path / 'index1'
    index2 = tmp_path / 'index2'
    ped1 = tmp_path / 'db1'
    ped2 = tmp_path / 'db2'
    index1.mkdir()
    index2.mkdir()
    ped1.mkdir()
    ped2.mkdir()
    # missing category column
    x = pd.DataFrame(
        {
            'giggle_index': [str(index1), str(index2)], 'ped_db': [str(ped1), str(ped2)]
        }
    )
    with pytest.raises(AssertionError):
        validate_stix_shardfile(x)
    # duplicates
    x = pd.DataFrame(
        {
            'giggle_index': [str(index1), str(index1)], 'ped_db': [str(ped1), str(ped1)], 'category': ['cat1', 'cat1']
        }
    )
    with pytest.raises(AssertionError):
        validate_stix_shardfile(x)
    # path does not exist
    x = pd.DataFrame(
        {
            'giggle_index': [str(tmp_path / 'missing_index')], 'ped_db': [str(ped1)], 'category': ['cat1']
        }
    )
    with pytest.raises(AssertionError):
        validate_stix_shardfile(x)

def test_read_giggle_shardfile_pass(tmp_path):
    index1 = tmp_path / 'index1'
    index2 = tmp_path / 'index2'
    index1.mkdir()
    index2.mkdir()

    shardfile_path = tmp_path / 'shardfile.tsv'
    with open(shardfile_path, 'w') as f:
        f.write('giggle_index\tcategory\n')
        f.write(f'{str(index1)}\tcat1\n')
        f.write(f'{str(index2)}\tcat2\n')

    df_giggle_shards = read_giggle_shardfile(str(shardfile_path), header=0)
    assert df_giggle_shards.shape[0] == 2
    assert set(df_giggle_shards['giggle_index']) == set([str(index1), str(index2)])
    assert set(df_giggle_shards['category']) == set(['cat1', 'cat2'])

def test_read_giggle_shardfile_fail(tmp_path):
    # index 2 path does not exist
    index1 = tmp_path / 'index1'
    index1.mkdir()

    shardfile_path = tmp_path / 'shardfile.tsv'
    with open(shardfile_path, 'w') as f:
        f.write('giggle_index\tcategory\n')
        f.write(f'{str(index1)}\tcat1\n')
        f.write(f'{str(tmp_path / "missing_index")}\tcat2\n')

    with pytest.raises(AssertionError):
        df_giggle_shards = read_giggle_shardfile(str(shardfile_path), header=0)


def test_read_stix_output_pass(tmp_path):
    # create a dummy stix output file
    stix_output_path = tmp_path / 'stix_output.tsv'
    with open(stix_output_path, 'w') as f:
        f.write("#gene_left=gene_x\n")
        f.write("#gene_right=gene_y\n")
        f.write("Miscellaneous header info\n")
        f.write('Giggle_File_Id\tPairend\tSplit\n')
        f.write('0\t1\t2\n')
        f.write('1\t2\t4\n')

    # this calls validate_stix_output, which will has an assertion check for columns
    gene_left, gene_right, df_stix_output = read_stix_fusion_output(str(stix_output_path))

def test_read_stix_output_fail(tmp_path):
    # create a dummy stix output file missing split column
    stix_output_path = tmp_path / 'stix_output.tsv'
    with open(stix_output_path, 'w') as f:
        f.write("#gene_left=gene_x\n")
        f.write("#gene_right=gene_y\n")
        f.write("Miscellaneous header info\n")
        f.write('Giggle_File_Id\tPairend\n')
        f.write('0\t1\n')
        f.write('1\t2\n')

    with pytest.raises(AssertionError):
        gene_left, gene_right, df_stix_output = read_stix_fusion_output(str(stix_output_path))
