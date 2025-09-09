## goal

Compute fraction of haas fusions that I have intervals for on GRCh37

### approach

- get set of haas genes
    - count num. occurences that each gene occurs in a fusion
        - gene_a: 1, gene_b: 3, ... 
- get set of genes from my latest bed file
- take the intersection
    - use haas gene count to get fraction of fusions that I have interval data for

### expected outcomes

- >80% good to run 1kg experiment
- >50% & <80% still good to run, but will need to add more intervals later
- <50% don't run the experiment yet


### code

```
# make maps of genes to aliases
# hgnc
./build_hgnc_alias_map.py
# haas
./build_haas_alias_map.py
```

- next step, write code to check if fusions are in the bed file I have, given the alias lookups

