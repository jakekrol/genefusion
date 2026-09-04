import pandas as pd
from importlib import resources

from importlib.resources import files

def get_data_file(filename):
    return files("polymerization.data").joinpath(filename)

def get_babiceanu_recurrent_normal_tissue_specific_fusions():
	path = get_data_file("babiceanu_recurrent_normal_tissue_specific_fusions.tsv")
	return pd.read_csv(path, sep="\t")

def get_babiceanu_recurrent_normal_tissue_agnostic_fusions():
	path = get_data_file("babiceanu_recurrent_normal_tissue_agnostic_fusions.tsv")
	return pd.read_csv(path, sep="\t")

def get_pcawg_recurrent_tumor_fusions():
	path = get_data_file("pcawg_recurrent_tumor_fusions.tsv")
	return pd.read_csv(path, sep="\t")

def get_pcawg_tumor_fusions():
	path = get_data_file("pcawg_fusions.tsv")
	return pd.read_csv(path, sep="\t")

def get_cosmic_tumor_fusions():
	path = get_data_file("cosmic_fusions.tsv")
	return pd.read_csv(path, sep="\t")

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
