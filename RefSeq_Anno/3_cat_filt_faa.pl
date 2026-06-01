use strict;
use warnings;

opendir DIR,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/faa_data" or die $!;
my @samples=readdir(DIR);
closedir DIR;

open OUTPUT,">/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/refseq_bacteria_assembly_summary_20240702_Salmonella_cds.faa0" or die $!;
foreach my $sample(@samples){
	next if($sample eq ".");
	next if($sample eq "..");
	my @arr=split(/_/,$sample);
	my $name=$arr[0]."_".$arr[1];
	open INPUT,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/faa_data/$sample" or die "cannot oepn faa_data/$sample\n";
	while(<INPUT>){
		if(/^>(.*)/){
			print OUTPUT ">$name|$1\n";
		}else{
			print OUTPUT $_;
		}
	}
	close INPUT;
}
close OUTPUT;
