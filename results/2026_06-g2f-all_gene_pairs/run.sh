#!/usr/bin/env bash

SEED=0
SAMPLE_SIZE=100000
GENOME_LIB_DIR=/data/jake/FusionAnnotator/genome_lib_dir
DIRBAM=../bams

# get fusion evidence
./run_g2f.py

# aggregate evidence across genes within modality
./run_agg.py

# stack modalitities
./join_modalities.sh

# split fusion evidence into chunks
./chunk.sh

# score
./score_chunks.py

# recollect
./combine_chunks.sh

# remove self-fusions
cd /data/jake/genefusion/results/2026_06-g2f-all_gene_pairs/chunks_scored
files=(blood_scored.tsv.sorted
bone_scored.tsv.sorted
breast_scored.tsv.sorted
esophagus_scored.tsv.sorted
gallbladder_scored.tsv.sorted
headneck_scored.tsv.sorted
kidney_scored.tsv.sorted
liver_scored.tsv.sorted
ovary_scored.tsv.sorted
pancreas_scored.tsv.sorted
prostate_scored.tsv.sorted)

for f in "${files[@]}"; do
    out="${f}.no_selfie"
    echo "# processing $f to outfile: $out"
    printf "gene_left\tgene_right\tscore_uniform\n" > $out
    tail -n +2 $f | awk '$1 != $2' >> $out
done

# plot tissue-wise score distributions
files=(blood_scored.tsv.sorted.no_selfie
bone_scored.tsv.sorted.no_selfie
breast_scored.tsv.sorted.no_selfie
esophagus_scored.tsv.sorted.no_selfie
gallbladder_scored.tsv.sorted.no_selfie
headneck_scored.tsv.sorted.no_selfie
kidney_scored.tsv.sorted.no_selfie
liver_scored.tsv.sorted.no_selfie
ovary_scored.tsv.sorted.no_selfie
pancreas_scored.tsv.sorted.no_selfie
prostate_scored.tsv.sorted.no_selfie)

for f in "${files[@]}"; do
	echo "# plotting $f"
	outdata="${f}.hist"
	outpng="${f}.hist.png"
	tissue=$(echo $f | cut -d'_' -f1)
	./hist_scores.py \
		--score "$f" \
		--seed "$SEED" \
		--sample_size "$SAMPLE_SIZE" \
		--out_data "$outdata" \
		--out_png "$outpng" \
		--tissue "$tissue"
done

# get top 10k fusions per tissue
for f in "${files[@]}"; do
	out="${f}.top10k"
	echo "# processing $f to outfile: $out"
	head -n 10001 $f > $out
done

# FusionAnnotator
files=(blood_scored.tsv.sorted.no_selfie.top10k
bone_scored.tsv.sorted.no_selfie.top10k
breast_scored.tsv.sorted.no_selfie.top10k
esophagus_scored.tsv.sorted.no_selfie.top10k
gallbladder_scored.tsv.sorted.no_selfie.top10k
headneck_scored.tsv.sorted.no_selfie.top10k
kidney_scored.tsv.sorted.no_selfie.top10k
liver_scored.tsv.sorted.no_selfie.top10k
ovary_scored.tsv.sorted.no_selfie.top10k
pancreas_scored.tsv.sorted.no_selfie.top10k
prostate_scored.tsv.sorted.no_selfie.top10k)
for f in "${files[@]}"; do
	out=${f}.annotated
	tmpfile=$(mktemp)
	trap 'rm -f "$tmpfile"' EXIT
	echo "# running FusionAnnotator on $f"
	printf "gene_left\tgene_right\tscore_uniform\n" > $tmpfile
	tail -n +2 $f | \
		awk -v OFS='\t' '{print $1"--"$2, $1, $2, $3}' >> $tmpfile
	FusionAnnotator \
		--genome_lib_dir $GENOME_LIB_DIR \
		--annotate $tmpfile | \
		grep -v PARA | grep -v NEIGHBOR | grep -v BLASTPAIR | grep -i oncogene > \
		$out
	rm -f "$tmpfile"
done

# subset to top 50
files=(blood_scored.tsv.sorted.no_selfie.top10k.annotated
bone_scored.tsv.sorted.no_selfie.top10k.annotated
breast_scored.tsv.sorted.no_selfie.top10k.annotated
esophagus_scored.tsv.sorted.no_selfie.top10k.annotated
gallbladder_scored.tsv.sorted.no_selfie.top10k.annotated
headneck_scored.tsv.sorted.no_selfie.top10k.annotated
kidney_scored.tsv.sorted.no_selfie.top10k.annotated
liver_scored.tsv.sorted.no_selfie.top10k.annotated
ovary_scored.tsv.sorted.no_selfie.top10k.annotated
pancreas_scored.tsv.sorted.no_selfie.top10k.annotated
prostate_scored.tsv.sorted.no_selfie.top10k.annotated)
for f in "${files[@]}"; do
	out=${f}.top50
	printf "fusion\tgene_left\tgene_right\tscore_uniform\tannotation\n" > $out
	echo "# processing $f to outfile: $out"
	head -n 50 $f >> $out
done

# get samples for each top fusion, for inspection
tissue_w_rna=(blood kidney liver ovary)
for tissue in "${tissue_w_rna[@]}"; do
	infile=${tissue}_scored.tsv.sorted.no_selfie.top10k.annotated.top50
	./top_fusions2samples.py \
		--input $infile \
		--tissue $tissue \
		--g2f_outdir ../g2f_out \
		--output ${infile}.samples
done

# download bams
for tissue in "${tissue_w_rna[@]}"; do
	infile=${tissue}_scored.tsv.sorted.no_selfie.top10k.annotated.top50.samples
	outfile=${infile}.bams.tsv
	./download_region_bams.py \
		--input $infile \
		--tissue $tissue \
		--output_dir $DIRBAM \
		--out_tbl $outfile \
		--bucket "layerlabcu"
done

# samplot

