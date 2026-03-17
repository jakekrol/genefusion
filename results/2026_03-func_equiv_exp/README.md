# Goal

Demonstrate gene connectivity is a reasonable approach to quantifying recurrent and functionally equivalent structural variants

## Data

Prostate cancer ERG-TMPRSS2 stix queries

## Analysis

How does # samples with ERG-TMPRSS2 fusion change as breakpoint merge thresholds are relaxed?

## Details

Sample-wise breakpoints are defined as (l,r). L is found by taking the minimum left read start coordinate among all read pairs supporting the fusion in the sample. R is found by taking maximum right read coordinate among all read pairs supporting the fusion in the sample.
Distance between sample i,j is the euclidean distance of (li,ri) and (lj,rj).
Compute all pairwise dij then threshold from dt=0 … d=max(dij) and for each dt count the number of samples with functional equiv SVs.