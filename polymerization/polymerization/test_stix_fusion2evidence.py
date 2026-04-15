import pytest
from polymerization.io import *
from polymerization.stix2fusion import *

def test_stix_fusion2evidence_pass(tmp_path):
    # create a dummy stix output file
    stix_output_path = tmp_path / 'stix_output.tsv'
    with open(stix_output_path, 'w') as f:
        f.write("#gene_left=gene_x\n")
        f.write("#gene_right=gene_y\n")
        f.write("Miscellaneous header info\n")
        f.write('Giggle_File_Id\tPairend\tSplit\n')
        f.write('0\t1\t2\n')
        f.write('1\t2\t4\n')

    gene_left, gene_right, df_stix_output = read_stix_fusion_output(stix_output_path)
    reads, samples = stix_fusion2evidence(df_stix_output)
    print(f"Reads: {reads}, Samples: {samples}")
    assert reads == 9
    assert samples == 2

def test_stix_fusion2evidence_fail(tmp_path):
	# negative evidence fails
	stix_output_path = tmp_path / 'stix_output.tsv'
	with open(stix_output_path, 'w') as f:
		f.write("#gene_left=gene_x\n")
		f.write("#gene_right=gene_y\n")
		f.write("Miscellaneous header info\n")
		f.write('Giggle_File_Id\tPairend\tSplit\n')
		f.write('0\t-5\t1\n')
		f.write('1\t1\t1\n')
	gene_left, gene_right, df_stix_output = read_stix_fusion_output(stix_output_path)

	with pytest.raises(AssertionError):
		reads, samples = stix_fusion2evidence(df_stix_output)


    