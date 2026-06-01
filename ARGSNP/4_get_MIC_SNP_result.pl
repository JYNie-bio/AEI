use strict;
use warnings;

my ($input1,$input2,$outfile)=@ARGV;

my %hash;
open INPUT1,"$input1" or die "cannot open $input1\n";
while(<INPUT1>){
	chomp;
	my @arr=split(/\t/);
	my $str=$arr[2];
	for(my $i=3;$i<=$#arr;$i++){
		$str=$str."\t".$arr[$i];
	}
	$hash{$arr[0]}=$str;
	if($arr[1] ne ""){
		$hash{$arr[1]}=$str;
	}
}
close INPUT1;

open OUTPUT,">$outfile" or die "cannot write $outfile\n";
open INPUT2,"$input2" or die "cannot open $input2\n";
my $line=<INPUT2>;
print OUTPUT "Genome_name\tCIP\tCHL\tNAL\tGEN\tTET\tCTX\tCFX\tAMP\tAMS\tCAZ\tCFZ\tIPM\tAZM\tSXT\tCT\t$line";
while(<INPUT2>){
	chomp;
	my @arr=split(/\t/);
	if(exists $hash{$arr[0]}){
		print OUTPUT $arr[0]."\t".$hash{$arr[0]}."\t".$_."\n";
	}
}
close INPUT2;
close OUTPUT;
