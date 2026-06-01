use strict;
use warnings;

my $temp="_add_mutation_models";
if($ARGV[0]=="1"){
	$temp="";
}
my %hash;
my %proteins;
open INPUT,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg${temp}_last_use_result_last.txt";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Sample");
	if(exists $hash{$arr[0]}{$arr[45]}){
		$hash{$arr[0]}{$arr[45]}++;
	}else{
		$hash{$arr[0]}{$arr[45]}=1;
	}
	$proteins{$arr[45]}=1;
}
close INPUT;

open OUTPUT,">/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/14880file_ARGanno_Protein_result_copynumber${temp}.txt";
open OUTPUT2,">/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/14880file_ARGanno_Protein_result_absence01${temp}.txt";
print OUTPUT "Genome";
print OUTPUT2 "Genome";
foreach my $prot(sort keys %proteins){
	print OUTPUT "\t".$prot;
	print OUTPUT2 "\t".$prot;
}
print OUTPUT "\n";
print OUTPUT2 "\n";
foreach my $sample(sort keys %hash){
	print OUTPUT $sample;
	print OUTPUT2 $sample;
	foreach my $prot2(sort keys %proteins){
		if(exists $hash{$sample}{$prot2}){
			print OUTPUT "\t".$hash{$sample}{$prot2};
			print OUTPUT2 "\t1";
		}else{
			print OUTPUT "\t0";
			print OUTPUT2 "\t0";
		}
	}
	print OUTPUT "\n";
	print OUTPUT2 "\n";
}
close OUTPUT;
close OUTPUT2;
