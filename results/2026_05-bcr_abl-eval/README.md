## goal:

show BCR-ABL would be prioritized if read depths were low, but recurrent in a population of tumor samples

## data

K562 and MCF-7 RNA-seq from CCLE and HAAS benchmark paper

## metric

Log10(p-valueBCR-ABL) under score distribution

## variation
- Samples
    - 1, …, 5 K562 samples
    - MCF-7 = 5 - K562 samples
- Reads
    - 20, 40, …, 100% sampled reads

## ideal result

Log10(p-valueBCR-ABL) high for low read depths, but high sample appearance.

## implementation steps

1. Data

- Tumor
- https://zenodo.org/records/13363154
- K562
    - SRR521460_1.fastq.20M.fq.gz
    - SRR521460_2.fastq.20M.fq.gz
- MCF-7
    - MCF7.Left.fq.gz
    - MCF7.Right.fq.gz
- Control/normal
    - 1000 genomes high coverage
- Sample K562 and MCF-7 at (20, 40, …, 100)

2. Run exhaustive g2f

- K562 (all samplings)
- MCF-7 (all samplings)
- 1000G (constant)

3. Score

- For each read-depth sampling
    - For each K562 sample fraction in {⅕, ⅖, …, 1}
        - Compute log10(pBCR-ABL) under score distribution