#!/bin/sh
outdir=$1;
reference=$2;
inputfile=$3;
faprefix=$4;
thread=$5;

date
cd ${outdir}
~/software/mummer-4.0.0rc1/nucmer --delta=${outdir}/${faprefix}.delta -t $thread $reference $inputfile
~/software/mummer-4.0.0rc1/show-snps -l -q -T ${faprefix}.delta >${faprefix}.snps0
~/software/mummer-4.0.0rc1/show-coords -rcl ${faprefix}.delta>${faprefix}.delta.coords
perl ~/ARGSNP/1_0_callsnp_filter_snps_nodeltafilter.pl ${faprefix}.snps0 ${faprefix}.delta.coords ${faprefix}.snps 95 90 90
perl ~/ARGSNP/1_1_get_mummcer_snp_indel_nodeltafilter.pl ~/MGE/MGE_result/Salmonella ${faprefix}.snps ${faprefix}.snps.snp.indel959090

date
