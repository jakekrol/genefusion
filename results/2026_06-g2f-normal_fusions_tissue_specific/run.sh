#!/usr/bin/env bash

# conda activate polymerization

### setup inputs
./get_g2f_inputs.sh

./sort_query_fusions.py

### get fusion evidence
./run_g2f.py

### aggregate
./join_evidence.sh

### get tissue specificity of normal fusions
./fusion2tissue.py

# ### score
# ./score_tumor_fusions.py

# ### plot
# mapfile -t files < <(ls tumor_score.* | grep -v png)
# for f in "${files[@]}"; do
#     tail -n +2 $f | \
#         bars.py -o ${f}.png --title $f
# done