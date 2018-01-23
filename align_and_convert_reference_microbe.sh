#!/bin/bash

echo $1 #first argument should be the base file name
echo "Generating GSNAPL output and phylowalked m8 files for sample:"
echo $1.fasta


echo "Running GSNAPL-1000"
/data/bin/gsnap-2015-12-31-bin/bin/gsnapl -A m8 --batch=4 --npaths=1000 --nofails --ordered -t 16 --maxsearch=1001 --max-mismatches=20 -D /nvme/databases/gsnap/nt-Jul2015-mask350-k16/nt_mask350-k16/ -d nt_mask350-k16 $1.fasta > $1.1000.m8
echo "Running phylowalk 1000"
python3 ./scripts/phylowalk_m8parser.py -i $1.1000.m8 -o $1.1000.m8.lca
