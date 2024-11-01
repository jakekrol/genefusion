#!/usr/bin/env bash

echo "begin"
date
# -B '/data/jake:/data' 
cmd="singularity shell --bind /data:/data ./stix-env.sif"
echo "running '$cmd'"
$cmd
echo "end"
date