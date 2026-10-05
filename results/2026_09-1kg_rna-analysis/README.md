1000G RNA MAGE short-read fusion analysis

1. What are the most common fusions called by star fusion? (bar top 10 and hist all)
2. How many fusions called per sample? What is the ratio of calls to chimeric reads? (hist and hist)
3. What are the highest scoring fusions with our approach? How does it compare to star fusion calls? (hist with star fusion calls as vlines)

discussion

1. Top 3 fusions were,

- TIMM23-PARGP1: published as a tumor progression fusion https://www.nature.com/articles/s41419-025-08106-w#additional-information
- LEPR-LEPROT: Leptin Receptor Overlapping Transcript. These sites are overlapping
- KLF16-OAZ1: reported as a hereditary fusion gene. 67.9% acute myeloid leukemia 	12.8% multiple myeloma 	5.2% GTEx blood (control)

another interesting top 10 recurrent fusion was: rom which KANSARL (KANSL1- ARL17A) was validated as the first predisposition fusion gene specific to 29% of populations of European ancestry.

2. Fusions called per sample. ratio of calls to chimeric reads.

- Fusions called per sample ususally less than 20
- Number of chimeric reads is pretty constant ~ 1M to 10M and has no relationship to num. calls

3. STIX is probably too sensitive with 1 read threshold per sample
