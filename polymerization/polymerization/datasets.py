import pandas as pd
from importlib import resources


def get_recurrent_normal_tissue_specific_fusions():
	"""
	Load recurrent normal tissue-specific fusions data as a pandas DataFrame.
	
	Returns:
		pd.DataFrame: DataFrame with columns: gene_left, gene_right, tissues
	"""
	# Load the TSV file from the package data directory
	if hasattr(resources, 'files'):
		# Python 3.9+
		data_file = resources.files('polymerization').joinpath('data', 'recurrent_normal_tissue_specific_fusions.tsv')
		with resources.as_file(data_file) as path:
			df = pd.read_csv(path, sep='\t')
	else:
		# Fallback for older Python versions
		data = resources.read_text('polymerization.data', 'recurrent_normal_tissue_specific_fusions.tsv')
		from io import StringIO
		df = pd.read_csv(StringIO(data), sep='\t')
	
	return df

def get_pcawg_recurrent_tumor_fusions():
	"""
	Load PCAWG recurrent tumor fusions data as a pandas DataFrame.
	
	Returns:
		pd.DataFrame: DataFrame with columns: gene_x, gene_y, projects, tissues
	"""
	# Load the TSV file from the package data directory
	if hasattr(resources, 'files'):
		# Python 3.9+
		data_file = resources.files('polymerization').joinpath('data', 'recurrent_tumor_fusions.tsv')
		with resources.as_file(data_file) as path:
			df = pd.read_csv(path, sep='\t')
	else:
		# Fallback for older Python versions
		data = resources.read_text('polymerization.data', 'recurrent_tumor_fusions.tsv')
		from io import StringIO
		df = pd.read_csv(StringIO(data), sep='\t')
	
	return df

def get_pcawg_data_types():
	x =	{
			'blood': 
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': ['tumor']
							}
						}
				},
			'bone': 
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': []
							}
						}
				},
			'breast': 
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': []
							}
						}
				},
			'esophagus':
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': []
							}
						}
				},
			'gallbladder':
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': []
							}
						}
				},
			'headneck':
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': []
							}
						}

				},
			'kidney':
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': ['tumor', 'normal']
							}
						}
				},
			'liver':
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': ['tumor', 'normal']
							}
						}
				},
			'ovary':
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': ['tumor']
							}
						}
				},
			'pancreas':
				{
					'modality':
						{
							'dna': {
								'specimen': []
							},
							'rna': {
								'specimen': ['tumor']
							}
						}
				},
			'prostate':
				{
					'modality':
						{
							'dna': {
								'specimen': ['tumor', 'normal']
							},
							'rna': {
								'specimen': []
							}
						}
				},
	}
	return x
