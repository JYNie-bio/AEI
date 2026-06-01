use strict;
use warnings;

my ($inputinfo,$inputlineage,$outfile)=(@ARGV);
#inputinfo: ~/Feature_analysis/updownstreamInfo/all_species_MGEresult_crossSpecies_prop.txt
#inputlineage: ~/MGE/bact.lineage2.csv
#outfile: ~/Feature_analysis/updownstreamInfo/updownstreamInfo/all_species_MGEresult_crossSpecies_prop_addlineage.txt

my %hash;
open INPUT0,"$inputlineage" or die $!;
while(<INPUT0>){
	chomp;
	my @arr=split(/\,/);
	$arr[8]=~s/ /_/g;
	$hash{$arr[8]}="$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[5]\t$arr[6]\t$arr[7]";	
}
close INPUT0;

open INPUT,"$inputinfo" or die $!;
open OUTPUT,">$outfile" or die $!;
my $line=<INPUT>;
chomp($line);
print OUTPUT "$line\tsuperkingdom\tkingdom\tphylum\tclass_bio\torder\tfamily\tgenus\n";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	if(!exists $hash{$arr[2]}){
		print "ERROR: donnot have species $arr[2] in $inputlineage\n";
	}else{
		print OUTPUT $_."\t".$hash{$arr[2]}."\n";
	}
}
close INPUT;
close OUTPUT;
