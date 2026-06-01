use strict;
use warnings;

my ($inputfile0,$inputfile1,$outfile)=(@ARGV);
#inputfile0: ~/MGE/all_species_runned_genome_ANI96need_genomelist.txt
#inputfile1: ~/Feature_analysis/updownstreamInfo/all_species_MGEresult_filtANI96_MobileElementFinder_BacAnt_PlasmidFinder.txt
#utfile: ~/Feature_analysis/updownstreamInfo/all_species_MGEresult_filtANI96_MobileElementFinder_BacAnt_PlasmidFinder_crossSpecies_prop.txt

my %hashnum;
open INPUT0,"$inputfile0" or die "cannot open $inputfile0\n";
while(<INPUT0>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Genome");
#	if(! exists $hashnum{$arr[0]}){
#		$hashnum{$arr[0]}=1;
#	}else{
#		$hashnum{$arr[0]}++;
#	}
	if(! exists $hashnum{$arr[2]}){
		$hashnum{$arr[2]}=1;
	}else{
		$hashnum{$arr[2]}++;
	}
}
close INPUT0;

my %hash;
open INPUT,"$inputfile1" or die "cannot open $inputfile1\n";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Species");
	my $str="$arr[3]\t$arr[10]";
	if(! exists $hash{$arr[0]}{$str}){
		$hash{$arr[0]}{$str}=$arr[1];
	}else{
		if($hash{$arr[0]}{$str}!~/\b\Q$arr[1]\E\b/){
			$hash{$arr[0]}{$str}=$hash{$arr[0]}{$str}.";".$arr[1];
		}
	}
}
close INPUT;

open OUTPUT,">$outfile" or die "cannot write $outfile\n";
print OUTPUT "MGE_Type\tElement_Type\tSpecies\tAnnoMGE_genome_num\tSpecies_allanno_num\tAnnoMGE_genome_num_proportion\n";
foreach my $species(sort keys %hash){
	foreach my $mge(sort keys %{$hash{$species}}){
		my @arr=split(/;/,$hash{$species}{$mge});
		my $annomgenum=$#arr+1;
		my $proportion=$annomgenum/$hashnum{$species};
		print OUTPUT $mge."\t".$species."\t".$annomgenum."\t".$hashnum{$species}."\t",$proportion."\n";
	}
}
close OUTPUT;
