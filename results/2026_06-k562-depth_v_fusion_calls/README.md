# goal

understand how depth affects number of fusion calls

# approach

- down/up-sample k562 cell line
- run star fusion on each depth
- count number of fusions called

# impact

- if some fusions are low depth, then single-sample callers may overlook them
- justifies the approach of aggregating across samples to call tumor fusions