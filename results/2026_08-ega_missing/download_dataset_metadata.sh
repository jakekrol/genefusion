#!/usr/bin/env bash
OUTDIR=ega_metadata
mkdir -p $OUTDIR
echo "# ega login"
egafetch-linux-amd64 auth login \
	--config-file ~/.egafetch/login.json
echo "# download ega dataset metadata"
for id in $(cat ega_datasets.txt); do
	echo "# $id"
	egafetch-linux-amd64 metadata $id \
		--format tsv \
		--config-file ~/.egafetch/login.json \
		--output ${OUTDIR}/${id}-metdata
done