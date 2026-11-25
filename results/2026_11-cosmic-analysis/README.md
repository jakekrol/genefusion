# goals

- Quantify COSMIC fusion evidence in 1KG and PCAWG normal
- understand expected fusion signal in non-tumor samples

## questions

- what is read/sample dist. for COSMIC fusions in 1KG & PCAWG normal, across DNA and RNA technologies
    - how does this compare to random gene pairs, not COSMIC.


## results

cosmic fusions with high evidence in 1kg and or pcawg normal

``` python
ANNOTATIONS_SCATTER=[
    ('GOPC', 'ROS1', 0.85, 0.15),
    ('NAB2', 'STAT6', 0.6, 0.9),
    ('CD74', 'NRG1', 0.2, 0.4)
]
```

GOPC--ROS1. intrachromosomal (chr6). overlapping. glioma (CNS) literature: GOPC:ROS1 and other ROS1 fusions represent a rare but recurrent drug target in a variety of glioma types. Acta Neuropathol. 2021 Dec;142(6):1065-1069. doi: 10.1007/s00401-021-02369-1. Epub 2021 Sep 18. PMID: 34536122; PMCID: PMC8568855.
NAB2--STAT6. intrachromosomal (chr12). adjacent. Solitary fibrous tumors:  NAB2–STAT6 fusion types account for clinicopathological variations in solitary fibrous tumors. Mod Pathol 28, 1324–1335 (2015). https://doi.org/10.1038/modpathol.2015.90
CD74--NRG1. interchromosomal (chr5, chr8). lung adenocarcinomas: CD74-NRG1 Fusions Are Oncogenic In Vivo and Induce Therapeutically Tractable ERBB2:ERBB3 Heterodimerization. Mol Cancer Ther 1 May 2022; 21 (5): 821–830. https://doi.org/10.1158/1535-7163.MCT-21-0820
