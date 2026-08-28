#!/usr/bin/env bash


CPUS=4
OUTDIR=bed
mkdir -p $OUTDIR

parallel_file="chimeric2bedpe.input"
mapfile -t chimeric_files < <(ls star_fusion-parallel/*/Chimeric.out.junction)
echo "${chimeric_files[@]}"
x=$(mktemp)
trap "rm $x" EXIT
for f in "${chimeric_files[@]}"; do
	echo $f >> $x
done

y=$(mktemp)
trap "rm $y" EXIT
sed "s|star_fusion-parallel/|$OUTDIR/|" $x > $y
sed -i "s|/Chimeric.out.junction|.bed.gz|" $y
paste $x $y > $parallel_file

chimeric2bedpe () {
	local chimeric_file="$1"
	local bedpe_file="$2"
	python -c "from polymerization.star import *; chimeric2bedpe(\"$chimeric_file\", \"$bedpe_file\", has_header=True, bgzip=True)"
}
export -f chimeric2bedpe

cat $parallel_file | gargs --log=chimeric2bedpe.log -p $CPUS "chimeric2bedpe {0} {1}"

