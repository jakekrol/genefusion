# goal

- for top scoring intrachromsomal fusions, check their supporting samples for read evidence

- workflow
    - select top fusions from score tables
    `get_top_candidates.sh`
    - re-query stix indices to get supporting samples
    - if num. supporting samples > 0
        - lookup associated rna file for sample
        - plot rna bam region
    
- note: only have rna seq for blood, esophagus, kidney, liver, and ovary.

Let's start with liver bc we have the most samples from there.