#!/bin/sh
faafile=$1
outdir=$2
faprefix=$3
#faprefix: GCA_010942655.1
date
cd ~/project/bjCDC/Feature_analysis
mkdir -p ${outdir}${faprefix}
sh 1_1_ARG_analyse5.sh $faafile ${outdir}${faprefix}/${faprefix}_ARG_combine_diamond.out ${outdir}${faprefix}/${faprefix}_ARG_combine_diamond.out_addpoint_models
sh 1_2_RefSeq_diamond.sh $faafile ${outdir}${faprefix}/${faprefix}_Salmonella_RefSeq_diamond.out
date
