#!/usr/bin/env bash

logfile='test_fmt_giggle_fusion.log'
echo "# test_fmt_giggle_fusion.sh"
echo "# logfile=$logfile"
cmd="gf-fmt_giggle_fusion -m test -g 7SK -r 8:144623280-144624570 -i data/7SK.8.pos.144623280.144624570"
echo "running '$cmd'"
$cmd > $logfile
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
cmd="gf-fmt_giggle_fusion -m test -g 7SK -r 8:144623280-144624570 -i data/7SK.8.pos.144623280.144624570 -o $tmp -z"
echo "running '$cmd'"
$cmd
echo "# test_fmt_giggle_fusion.sh done"