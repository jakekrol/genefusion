# goal

evaluate scores ability to classify tumor V normal fusions

# approach

for three different scoring techniques, use score as threshold in ROC 

i) aggregate evidence from all tissues, weighted by subpopulation size
ii) aggregate evidence from all tissues, uniform subpop weights
iii) score using only tissue where fusion occurs, uniform subpop-tissue weights