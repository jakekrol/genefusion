#!/usr/bin/env bash


CPUS=16
INDIR=bed
OUTDIR=bed_clean
mkdir -p $OUTDIR

parallel_file="clean_bedpe.input"
mapfile -t files_to_clean < <(ls ${INDIR}/*.bed.gz)
echo "${files_to_clean[@]}"
x=$(mktemp)
trap "rm $x" EXIT
for f in "${files_to_clean[@]}"; do
	echo $f >> $x
done

y=$(mktemp)
trap "rm $y" EXIT
sed "s|$INDIR/|$OUTDIR/|" $x > $y
paste $x $y > $parallel_file

# chimeric output bed is compatible with excord output format
clean_excord () {
	local infile="$1"
	local outfile="$2"
	python -c "from polymerization.giggle2fusion import *; clean_excord(\"$infile\", \"$outfile\", bgzip=True)"
}
export -f clean_excord

cat $parallel_file | gargs --log=clean_bedpe.log -p $CPUS "clean_excord {0} {1}"
