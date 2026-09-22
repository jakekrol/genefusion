grep BAM ../../data/2024_08-icgc_legacy_locations/icgc25k-legacy-data-locations.no_index.tsv \
    > icgc25-legacy-data.bam.tsv
./find_bams.py
./name_bams.py