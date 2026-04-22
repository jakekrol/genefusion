# goal

evaluate score approach

## data

- tissue-specific recurrent fusions in pcawg tumor samples
	- ../2026_04-s2f_pcawg_recurrent/all_categories-fusion_table.tsv
- tissue-agnostic recurrent fusions in normal samples
	- ../2026_04-s2f_babiceanu_recurrent_normal_tissue_agnostic/all_categories-fusion_table.tsv

## approach

- positives are always the recurrent tumor fusions
	- evidence types used for scoring
		- tumor (DNA and RNA when available) from recurrent tissue(s)
		- normal (DNA and RNA when available) from recurrent tissue(s)
		- optionally, high coverage 1000G
- negatives: 3 different negative sets
- re-use the recurrent tumor fusions, except score in tissues where they are not recurrent (off-tissue)
	- 1) score using all off-tissues (# negatives = # positives)
	- 2) score using a single off-tissue (# negatives = (# positives * # off-tissues)
- 3) negatives are recurrent normal fusions and score is computed using all tissues in PCAWG
- for each negative set, toggle 1000G inclusion/exclusion in score
- for each negative set experiment, vary the w_normal parameter {1/4, 1/2, 3/4}



		

1. positives are tissue-specific fusions. negatives are off-tissue agnostic fusions
2. positives are tissue-specific fusions. negatives are off-tissue specific fusions
3. positives are tissue-specific fusions. negatives are recurrent normal fusions using same evidence (tissue-specific) as positives

## Results

1. Score negatives are tissue agnostic (N=P)