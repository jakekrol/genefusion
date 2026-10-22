#!/usr/bin/env bash

CHUNKS=8
# fids=(EGAF00001719777 EGAF00001719891 EGAF00001719867)
fids=(EGAF00001719777 EGAF00001719867)
for f in "${fids[@]}"; do
    # continue until md5 exists. md5 is only written if download checksum is correct.
    while ! compgen -G "$f/*.md5" > /dev/null; do
        egafetch-linux-amd64 download $f \
            --config-file credentials.json --restart --parallel-chunks $CHUNKS
        sleep 3
    done
done
