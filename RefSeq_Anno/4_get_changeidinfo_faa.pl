use strict;
use warnings;

open INPUT,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/refseq_bacteria_assembly_summary_20240702_Salmonella_cds.faa0" or die $!;
open OUTPUT,">/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/refseq_bacteria_assembly_summary_20240702_Salmonella_cds.faa_changeid.info.txt" or die $!;
my $a=0;
while(<INPUT>){
	if(/^>/){
		$a++;
		print OUTPUT "id".$a."\t".$_;
	}
}
close INPUT;
close OUTPUT;
