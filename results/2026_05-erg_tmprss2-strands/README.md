# motivation

high read depth occurs in high cov 1000g for ERG-TMPRSS2 ~300 and BCR-ABL ~100. this is noise, and we wonder if strand configurations are biased in either tumor or normal. if yes, we can use this to filter/prioritize reads

# goal

- quantify strand configurations of ERG-TMPRSS2 fusion in tumor, paired normal, and 1000g

# approach

- run g2f for erg-tmprss2 on the 3 indices (prostaste tumor, prostate normal, and 1000g)
- scatter the supporting reads by breakpoint and color/fill by strand config.
- barplot strand config
