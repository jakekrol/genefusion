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
