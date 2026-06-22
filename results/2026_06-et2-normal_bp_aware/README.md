# goal

test out breakpoint aware normal read/sample counting in ERG-TMPRSS2

# approach

- bin 2d paired gene space
- only count normal reads that cancel out a neighboring tumor read (same bin)
	- mathematically $N_{ij}' = min(T_{ij},N_{ij})$

