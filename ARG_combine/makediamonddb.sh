#!/bin/sh
cd ~/ARG_combine
date
. ~/miniconda3/bin/activate diamond
cut -f1 ARG_combine_info_aa_use2.txt |sed -n '2,$p' |seqkit grep -f - ARG_combine_rmduplicated_addpoint_models.faa -o ARG_combine_rmduplicated.faa
diamond makedb --in ARG_combine_rmduplicated.faa -d ARG_combine_rmduplicated.faa.diamond
diamond makedb --in ARG_combine_rmduplicated_addpoint_models.faa -d ARG_combine_rmduplicated_addpoint_models.faa.diamond
date
