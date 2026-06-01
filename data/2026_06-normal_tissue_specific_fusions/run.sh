#!/usr/bin/env bash

url="https://docs.google.com/spreadsheets/d/1dzicCSrs1n17duJuIK8EGy1iSvLQ8JhW/export?format=tsv&gid=130561229"
curl -L -o recurrent_normal_tissue_specific.tsv "$url"

tail -n +2 recurrent_normal_tissue_specific.tsv > z && mv z recurrent_normal_tissue_specific.tsv