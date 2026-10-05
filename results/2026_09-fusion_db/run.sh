# documenting order of commands

# organize tables in sqlite
./create_db.py
# compute various scores for fusions
./score_chunks.py
# sort by score and select top fusions
./rank.py
# annotate and filter
./annotate.sh
# filter high repeat coverage
./repeat_filter.py
# estimate breakpoint
./estimate_breakpoint.py
# rank by breakpoint score
./rank_breakpoints.py
# find samples with maximal read evidence for top fusions
./topfusion2samples.py
# download sample bams for visual validation
./download_bams.py