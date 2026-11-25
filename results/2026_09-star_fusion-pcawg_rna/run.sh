grep BAM ../../data/2024_08-icgc_legacy_locations/icgc25k-legacy-data-locations.no_index.tsv \
    > icgc25-legacy-data.bam.tsv
./find_bams.py
./name_bams.py
./pipeline.py
./chimeric2bedpe.sh
./get_chimeric_counts_per_sample.sh
./get_star_fusion_calls.py
./clean_bedpe.sh
./sort_bed bed_clean bed_clean_sort
./sort_beds_by_tissue.py
