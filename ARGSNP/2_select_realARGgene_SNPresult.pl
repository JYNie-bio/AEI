use strict;
use warnings;

my ($inputfile,$samplename,$outfile)=(@ARGV);

my %hash;
my $inputfile0;
if($inputfile~/_add_mutation_models/){
	$inputfile0="~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last.txt";
}else{
	$inputfile0="~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_add_mutation_models_last_use_result_last.txt";
}
open INPUT0,$inputfile0 or die $!;
while(<INPUT0>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Sample");
	$hash{$arr[0]}{$arr[1]}="$arr[41]|$arr[44]|$arr[45]";
}
close INPUT0;

open INPUT,"$inputfile" or die "cannot open $inputfile\n";
open OUTPUT,">$outfile" or die "cannot open $outfile\n";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	if($arr[0] eq "Ref_contig"){
		print OUTPUT $_."\n";
	}else{
		if(exists $hash{$samplename}{$arr[4]}){
			my @array1=split(/\|/,$hash{$samplename}{$arr[4]});
			my @array2=split(/\|/,$arr[0]);
			if(($array1[1] eq $array2[3]) && ($array1[2] eq $array2[4])){
				my $str=$samplename."|".$arr[4]."|".$hash{$samplename}{$arr[4]};
				print OUTPUT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$str\t$arr[5]\t$arr[6]\t$arr[7]\t$arr[8]\t$arr[9]\t$arr[10]\t$arr[11]\t$arr[12]\t$arr[13]\t$arr[14]\n";
			}
		}
	}
}
close INPUT;
close OUTPUT;
