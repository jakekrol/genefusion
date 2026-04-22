#!/usr/bin/env bash


echo "# DEPRECATED. see plot_rocs-neg_tissue_agnostic.py"
exit 1
ROC_SCRIPT=/data/jake/rl-tools/plot/roc.py
scores=score.tsv

[[ -f $ROC_SCRIPT ]] || { echo "ROC script not found at $ROC_SCRIPT"; exit 1; }

# negative set 1
# negatives are recurrent tumor fusions, except scored in all tissues where they are not recurrent
# (off-tissue scoring)
w_normals=(0.0 0.25 0.5 0.75 1.0)
for w_normal in "${w_normals[@]}"; do
	outfile="roc_input-neg1-wnormal_$w_normal.tsv"
	script="
import pandas as pd
df=pd.read_csv('$scores', sep='\t')
cols_keep=['gene_left', 'gene_right']
p=df[cols_keep + ['score_as_positive-wnormal_$w_normal']]
p.columns = cols_keep + ['score']
p['label'] = 1
n=df[cols_keep + ['score_as_negative-wnormal_$w_normal']]
n.columns = cols_keep + ['score']
n['label'] = 0
df_out = pd.concat([p,n], axis=0)
df_out.sort_values('score', inplace=True, ascending=False)
df_out.to_csv('$outfile', sep='\t', index=False)
"
	python -c "$script"
	outfile_w_1000g="roc_input-neg1-wnormal_${w_normal}-w_1000g.tsv"
	script="import pandas as pd
df=pd.read_csv('$scores', sep='\t')
cols_keep=['gene_left', 'gene_right']
p=df[cols_keep + ['score_as_positive-w_1000g-wnormal_${w_normal}']]
p.columns = cols_keep + ['score']
p['label'] = 1
n=df[cols_keep + ['score_as_negative-w_1000g-wnormal_${w_normal}']]
n.columns = cols_keep + ['score']
n['label'] = 0
df_out = pd.concat([p,n], axis=0)
df_out.sort_values('score', inplace=True, ascending=False)
df_out.to_csv('$outfile_w_1000g', sep='\t', index=False)
"
	python -c "$script"	
done

mapfile -t roc_files < <(ls | grep roc | grep tsv)
roc_files="$(echo ${roc_files[@]} | sed 's/ /,/g')"
echo "ROC files: $roc_files"
$ROC_SCRIPT \
	--scores $roc_files \
	--output roc_curves.png \
	--header \
	--score_col 2 \
	--label_col 3 \
	--names "Wnormal=0,Wnormal=0.25,Wnormal=0.5,Wnormal=0.75,Wnormal=1" \
	--title "Recurrent tumor fusion classification"
