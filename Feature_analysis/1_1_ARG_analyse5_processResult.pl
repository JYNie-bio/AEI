#!/usr/bin/perl -w
use strict;
use Math::Round;
my ($faprefix)=@ARGV;
#$faprefix: ah19S225

my %hash;
open INPUT,"/lustre/home/niejingyi2023/database/ARG_combine/ARG_combine_info_aa.txt" or die $!;
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	$hash{$arr[0]}=$_;
}
close INPUT;

######## diamond_ARG_combine ######
open INPUT1,"/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/ARGcombine_RefSeq_result/${faprefix}/${faprefix}_ARG_combine_diamond.out" or die "Cannot open /lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/ARGcombine_RefSeq_result/${faprefix}/${faprefix}_ARG_combine_diamond.out\n";
open OUTPUT1,">/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/ARGcombine_RefSeq_result/${faprefix}/${faprefix}_ARG_combine_diamond.out.processed.txt" or die "Cannot write /lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/ARGcombine_RefSeq_result/${faprefix}/${faprefix}_ARG_combine_diamond.out.processed.txt\n";
print OUTPUT1 "query_id\tquery_length\tsubject_id\tsubject_length\talignment_length\tmismatches\tgap_opens\tqstart\tqend\tsstart\tsend\tidentity\tqcov\tscov\tevalue\tbit_score\tSubject_ID\tGene_name\tAlternative_name\tProduct_name\tGene_family\tDrug_class\tDrug_name\tMechanism_of_resistance\tOther_Notes\n";
while(<INPUT1>){
	chomp;
	next if(/^#/);
        my @arr=split(/\t/);
	my $qcov=100*(abs($arr[10]-$arr[9])+1)/$arr[1];
	my $scov=100*(abs($arr[12]-$arr[11])+1)/$arr[3];
	if(! exists $hash{$arr[2]}){
		print "ERROR: do not have $arr[2] in ARG_combine_info_aa.txt\n";
	}else{
		print OUTPUT1 "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[5]\t$arr[7]\t$arr[8]\t$arr[9]\t$arr[10]\t$arr[11]\t$arr[12]\t$arr[4]\t$qcov\t$scov\t$arr[13]\t$arr[14]\t$hash{$arr[2]}\n";
	}
}
close INPUT1;
close OUTPUT1;

