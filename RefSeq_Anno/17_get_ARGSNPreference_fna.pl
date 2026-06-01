use strict;
use warnings;

my ($inputfile,$outfile)=(@ARGV);
my %hash;
open INPUT,"$inputfile";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Sample");
	$hash{$arr[0]}{$arr[1]}="$arr[0]|$arr[1]|$arr[17]|$arr[18]|$arr[19]|$arr[8]|$arr[9]|$arr[10]|$arr[11]|$arr[12]|$arr[13]|$arr[14]|$arr[15]";
	#genome|protenid|Annotation_Source|Use_Gene_Name_Modify|Use_Protein_Name_Modify|alignment_length|mismatches|gap_opens|identity|qcov|scov|evalue|bit_score
}
close INPUT;

open OUTPUT,">$outfile" or die "cannot write\n";
foreach my $genome(sort keys %hash){
	my %hashffn=();
	my $str;
	open INPUT2,"/lustre/home/niejingyi2023/project/bjCDC/MGE/MGE_result/Salmonella/${genome}/prokka/${genome}.ffn" or die "cannot open /lustre/home/niejingyi2023/project/bjCDC/MGE/MGE_result/Salmonella/${genome}/prokka/${genome}.ffn\n";
	while(<INPUT2>){
		if(/^>(.*?) /){
			$str=$1;
			$hashffn{$str}="";
		}else{
			$hashffn{$str}=$hashffn{$str}.$_;
		}
	}
	close INPUT2;
	foreach my $protid(sort keys %hashffn){
		if(exists $hash{$genome}{$protid}){
			print OUTPUT ">".$hash{$genome}{$protid}."\n".$hashffn{$protid};
		}
	}
}
close OUTPUT;
