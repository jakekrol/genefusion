#!/usr/bin/env bash
# note: for some reason singularity cachedir cannot be set to data partition

#build with singularity build --sandbox --fakeroot stix-env.sif stix-env.def
echo "begin"
date
build_cmd="singularity build --sandbox --fakeroot stix-env.sif stix-env.def"
echo "running '$build_cmd'"
$build_cmd
date
echo "end"

#singularity shell -B '/data:/data' ~/Repositories/STIX_LR_paper/SR/stix-env.sif
#singularity shell -B '/data:/data/jake/sv' ~/Repositories/STIX_LR_paper/SR/stix-env.sif
