#!/bin/sh
faprefix=$1;
inputfile=~/MGE/MGE_result/Salmonella/${faprefix}/${faprefix}_MGE_merged_site.txt
awk -F'\t' '{if($3=="plasmid"){print $4}}' $inputfile>~/MGE/MGE_result/Salmonella_MOBsuite/ContigID/${faprefix}_plasmidcontigID.txt
seqkit grep -f ~/MGE/MGE_result/Salmonella_MOBsuite/ContigID/${faprefix}_plasmidcontigID.txt ~/MGE/MGE_result/Salmonella/${faprefix}/temp/${faprefix}.fa --threads 4 >~/MGE/MGE_result/Salmonella_MOBsuite/plamid_fasta/${faprefix}.plasmid.fa

