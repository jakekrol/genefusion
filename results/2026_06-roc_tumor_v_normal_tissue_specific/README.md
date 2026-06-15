# goal

evaluate scores ability to classify tumor V normal fusions

# approach

for three different scoring techniques, use score as threshold in ROC 

i) aggregate evidence from all tissues, weighted by subpopulation size
ii) aggregate evidence from all tissues, uniform subpop weights
iii) score using only tissue where fusion occurs, uniform subpop-tissue weights

# misclassification analysis

- ZNF765--TPM3P9
	- label: normal
	- score: ~0.2
	- data: {'reads_high_coverage_1000g_dna': 21.0, 'samples_high_coverage_1000g_dna': 13.0, 'reads_ovary_normal_dna': 32.0, 'samples_ovary_normal_dna': 17.0, 'reads_ovary_tumor_dna': 31.0, 'samples_ovary_tumor_dna': 21.0, 'reads_ovary_tumor_rna': 1946.0, 'samples_ovary_tumor_rna': 70.0}
	- reason: RNA-seq is only available for tumor samples, not paired normal, causing positive score
- SIDT2--TAGLN
	- label: normal
	- score: ~0.14
	- data:  {'reads_bone_normal_dna': 0.0, 'samples_bone_normal_dna': 0.0, 'reads_bone_tumor_dna': 8.0, 'samples_bone_tumor_dna': 2.0, 'reads_esophagus_normal_dna': 0.0, 'samples_esophagus_normal_dna': 0.0, 'reads_esophagus_normal_dna': 0.0, 'samples_esophagus_normal_dna': 0.0, 'reads_esophagus_tumor_dna': 0.0, 'samples_esophagus_tumor_dna': 0.0, 'reads_gallbladder_normal_dna': 0.0, 'samples_gallbladder_normal_dna': 0.0, 'reads_gallbladder_tumor_dna': 0.0, 'samples_gallbladder_tumor_dna': 0.0, 'reads_high_coverage_1000g_dna': 23.0, 'samples_high_coverage_1000g_dna': 14.0, 'reads_pancreas_tumor_rna': 691.0, 'samples_pancreas_tumor_rna': 68.0}
	- reason: RNA-seq is only available for tumor samples, not paired normal, causing positive score
- NDC80--SMCHD1
	- label: tumor
	- score: ~-0.006
	- data: {'reads_bone_normal_dna': 5.0, 'samples_bone_normal_dna': 4.0, 'reads_bone_tumor_dna': 3.0, 'samples_bone_tumor_dna': 1.0, 'reads_high_coverage_1000g_dna': 65.0, 'samples_high_coverage_1000g_dna': 30.0}
	- reason: low tumor depth, some depth in normal and 1000g
- STAT3--PTRF
	- label: tumor
	- score: ~-0.003
	- data:  {'reads_breast_normal_dna': 6.0, 'samples_breast_normal_dna': 3.0, 'reads_breast_tumor_dna': 1.0, 'samples_breast_tumor_dna': 1.0, 'reads_high_coverage_1000g_dna': 10.0, 'samples_high_coverage_1000g_dna': 5.0}
	- reason: low tumor depth, some depth in normal and 1000g

addressing problems

i) exclude RNA-seq
ii) see if breakpoint or strand config is diff for the low scoring tumor fusions


