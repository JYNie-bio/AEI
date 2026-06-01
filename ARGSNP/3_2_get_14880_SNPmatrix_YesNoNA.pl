use strict;
use warnings;

my ($snpfile,$mummerfolder,$outfile)=(@ARGV);
my %snps;
open INPUT,"$snpfile";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Sample");
	$snps{$arr[$#arr]}=1;
}
close INPUT;

opendir DIR,"$mummerfolder";
my @files=readdir(DIR);
closedir DIR;

my %hash;
my %samples;
my %hashgenome;
foreach my $file(@files){
	next if($file!~/snps.snp.indel959090.filtedARG$/);
	my $sample=$file=~s/.snps.snp.indel959090.filtedARG//r;
	$samples{$sample}=1;
	open INPUT2,"$mummerfolder/$file" or die "cannot open $mummerfolder/$file\n";
	while(<INPUT2>){
		chomp;
		my @arr=split(/\t/);
		next if($arr[0] eq "Ref_contig");
		my @arr2=split(/\|/,$arr[0]);
		my $str=$arr2[4]."|".$arr[3]."-".$arr[1]."-".$arr[2];
		$hash{$str}{$sample}=1;
		$hashgenome{$arr2[4]}{$sample}=1;
	}
	close INPUT2;
}

open OUTPUT,">$outfile";
print OUTPUT "Sample";
for my $key(sort keys %hash){
	print OUTPUT "\t$key";
}
print OUTPUT "\n";
for my $sample(sort keys %samples){
	print OUTPUT $sample;
	for my $key(sort keys %hash){
		my @array=split(/\|/,$key);
		if(exists $hash{$key}{$sample}){
			print OUTPUT "\tYes";
		}elsif(exists $hashgenome{$array[0]}{$sample}){
			print OUTPUT "\tNo";
		}else{
			print OUTPUT "\tNA";
		}
	}
	print OUTPUT "\n";
}
