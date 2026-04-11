import pytest
import pandas as pd
from polymerization.stix2fusion import *


def test_verify_fusion_set_in_bed_fail():
    df_fusion = pd.DataFrame(
        {
            'x': ['gene_a', 'gene_b'], 'y': ['gene_c', 'gene_d']
        }
    )
    df_bed = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr2'], 'start': [100, 150], 'end': [200, 300], 'gene_name': ['gene_y', 'gene_z']
        }
    )
    with pytest.raises(AssertionError):
        verify_fusion_set_in_bed(df_fusion, df_bed)

def test_verify_fusion_set_in_bed_pass():
    df_fusion = pd.DataFrame(
        {
            'x': ['gene_a', 'gene_b'], 'y': ['gene_c', 'gene_d']
        }
    )
    df_bed = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr2', 'chr3', 'chr4'], 'start': [100, 150, 200, 250], 'end': [200, 300, 400, 450], 'gene_name': ['gene_a', 'gene_b', 'gene_c', 'gene_d']
        }
    )
    verify_fusion_set_in_bed(df_fusion, df_bed)

def test_left_sort_fusion_set_pass():
    # test chromosome order
    df_fusion = pd.DataFrame(
        {
            'gene_x': ['gene_a'], 'gene_y': ['gene_b']
        }
    )
    df_bed = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr2'], 'start': [100, 50], 'end': [200, 100], 'gene_name': ['gene_a', 'gene_b']
        }
    )
    df_fusion_sorted = left_sort_fusion_set(df_fusion, df_bed)
    assert df_fusion_sorted.iloc[0]['gene_left'] == 'gene_a'
    assert df_fusion_sorted.iloc[0]['gene_right'] == 'gene_b'
    # test start position order
    df_fusion = pd.DataFrame(
        {
            'gene_x': ['gene_a'], 'gene_y': ['gene_b']
        }
    )
    df_bed = pd.DataFrame(
        {
            'chromosome': ['chr1', 'chr1'], 'start': [50, 100], 'end': [200, 200], 'gene_name': ['gene_a', 'gene_b']
        }
    )
    df_fusion_sorted = left_sort_fusion_set(df_fusion, df_bed)
    assert df_fusion_sorted.iloc[0]['gene_left'] == 'gene_a'
    assert df_fusion_sorted.iloc[0]['gene_right'] == 'gene_b'
