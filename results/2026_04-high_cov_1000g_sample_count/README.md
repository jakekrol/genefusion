count num samples in 1000g high cov index

```
# use polymerization environment for giggle
index_parent=/data/jake/1kg_high_coverage_2025
cd $index_parent
mapfile -t shards < <(ls | grep shard_ | grep -v ped)
for s in "${shards[@]}"; do
	echo $s
	giggle search -i $s -l | tail -n +2 | wc -l
done
# shard_0
# 313
# shard_1
# 313
# shard_2
# 313
# shard_3
# 313
# shard_4
# 313
# shard_5
# 313
# shard_6
# 313
# shard_7
# 313

```

313 * 8 is 2504