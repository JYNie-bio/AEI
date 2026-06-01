#!/bin/sh
faprefix=$1
threads=$2
date
##########The default thresholds are heavily optimized for publicly available Enterobacteriaceae plasmids and these may not be appropriate for other taxa of interest.

inputfa=~/MGE/MGE_result/Salmonella_MOBsuite/plamid_fasta/${faprefix}.plasmid.fa
outdir=~/MGE/MGE_result/Salmonella_MOBsuite/MOBsuite_result
outfile="${outdir}/${faprefix}_mobsuite.txt"
mob_typer -i $inputfa -o $outfile -n $threads -s $faprefix --multi
date
