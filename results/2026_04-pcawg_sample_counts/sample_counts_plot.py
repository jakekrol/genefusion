#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt

f = 'sample_counts.tsv'
data = pd.read_csv(f, sep='\t')
data['tissue'] = data['tissue'].str.capitalize()
data['modality'] = data['modality'].str.upper()

# Define combo labels and mellow colors
order = ['dna-tumor', 'dna-normal', 'rna-tumor', 'rna-normal']
colors = {
	'DNA-tumor': 'blue',
	'DNA-normal': 'lightblue',
	'RNA-tumor': 'red',
	'RNA-normal': "lightcoral",
}

data['group'] = data['modality'] + '-' + data['specimen']

# First plot data: keep modality as x-axis and specimen grouped within each modality
pivot_modality = data.pivot_table(index='modality', columns='specimen', values='samplecount', aggfunc='sum', fill_value=0)
pivot_modality = pivot_modality.reindex(index=['DNA', 'RNA'], fill_value=0)

# Create the first bar plot
plt.figure(figsize=(10, 5))

# Bar plot with tumor/normal adjacent within each modality
ax = plt.gca()
x = range(len(pivot_modality.index))
bar_width = 0.35

for i, modality in enumerate(pivot_modality.index):
	ax.bar(
		i - bar_width / 2,
		pivot_modality.loc[modality, 'normal'],
		width=bar_width,
		color='black',
		label=f'{modality.upper()}-normal',
	)
	ax.bar(
		i + bar_width / 2,
		pivot_modality.loc[modality, 'tumor'],
		width=bar_width,
		color='gray',
		label=f'{modality.upper()}-tumor',
	)

plt.title('PCAWG sample counts')
plt.ylabel('Sample count')
plt.xlabel('Modality')
plt.xticks(list(x), list(pivot_modality.index), rotation=0)
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), title='Group')
plt.tight_layout()
plt.savefig('sample_counts_modality.png')
plt.show()

# 2. Plot stratified by tissue
# Create a pivot table for tissue stratification (DNA + RNA)
pivot_tissue = data.pivot_table(index='tissue', columns='group', values='samplecount', aggfunc='sum', fill_value=0)
# Create the second bar plot
plt.figure(figsize=(12, 6))

# Bar plot for tissue stratification
pivot_tissue.plot(kind='bar', width=0.85, color=[colors[col] for col in pivot_tissue.columns])
plt.title('PCAWG sample counts stratified by tissue')
plt.ylabel('Sample count')
plt.xlabel('Tissue')
plt.xticks(rotation=45)
plt.legend(title='Modality')
plt.tight_layout()
plt.savefig('sample_counts_tissue.png')
plt.show()