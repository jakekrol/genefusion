#!/usr/bin/env bash
# parallel stix queries for pcawg prostate gene fusion data

procs=${1:-5}

gargs --log genefusion_stix_sharded.log -p $procs -o "genefusion-genefusion_stix_sharded {0} {1} {2} {3} {4}" < genefusion-stix_sharded.input