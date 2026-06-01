use strict;
use warnings;

open INPUT2,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/refseq_bacteria_assembly_summary_20240702_Salmonella_cds.faa_changeid.info.txt" or die $!;
my %hash;
while(<INPUT2>){
	chomp;
	my @arr=split(/\t/);
	$hash{$arr[1]}=$arr[0];
}
close INPUT2;

open INPUT,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/refseq_bacteria_assembly_summary_20240702_Salmonella_cds.faa0" or die $!;
open OUTPUT,">/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/refseq_bacteria_assembly_summary_20240702_Salmonella_cds.faa" or die $!;
while(<INPUT>){
	chomp;
	if(/^>/){
		if(exists $hash{$_}){
			print OUTPUT ">$hash{$_}\n";
		}else{
			print "ERROR: do not have $_ in faa0 file\n";
		}
	}else{
		print OUTPUT $_."\n";
	}
}
close OUTPUT;
close INPUT;
