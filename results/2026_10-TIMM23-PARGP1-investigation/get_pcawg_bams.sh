#!/usr/bin/env bash

CHUNKS=8
fids=(EGAF00001719777 EGAF00001719891 EGAF00001719867)
for f in "${fids[@]}"; do
    egafetch-linux-amd64 download $f \
        --config-file credentials.json --restart --parallel-chunks $CHUNKS
done
