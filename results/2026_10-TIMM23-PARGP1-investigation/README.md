# goal

compare TIMM23-PARGP1 in tumor and normal RNA-seq

## approach

- find 3 tumor and 3 normal RNA-seq samples with high read evidence
- samplot
- breakpoint plot
- find most likely transcripts in tumor and normal

## expected result

case i) tumor and normal are similar: breakpoint, transcript, and samplot

case ii) tumor and normal are dissimlar: breakpoint, transcript, and samplot


## misc.

PARGP1 transcript from GENCODE is near 50kbp. but my bed annotation is near 46kbp.

This inconsitency occurs in many DBs:

- GeneCard: ~46kbp
- Ensemble: ~50kbp
- NCBI: ~46kbp

I used th hg19.ncbiRefSeq.gtf.gz   GTF from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/
to get exons near 46kbp

mapping pcawg file ids to ega ids using pcawg_fileid2ega_fileid.tsv

| pcawg | ega |
| ----- | ----- |
| FI6523 | EGAF00001719777 |
| FI1218 | EGAF00001719891 |
| FI40382 | EGAF00001719867 |

