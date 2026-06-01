use strict;
use warnings;

my $filename=$ARGV[0];
my %infohash;
my %wphash;
open INPUT00,"~/RefSeq_Anno/RefSeq_Anno_Result_14888file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_needChangeInfo_add_mutation_models.txt";
while(<INPUT00>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Annotation_Source");
	my $str=$arr[0]."\t".$arr[1]."\t".$arr[2];
	$infohash{$str}=$arr[3]."\t".$arr[4];
	if($arr[3] eq "NA"){
		$wphash{$arr[2]}=1;
	}
}
close INPUT00;


open INPUT,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/${filename}_last_use_result.txt";
open OUTPUT,">/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/${filename}_last_use_result_last.txt";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	if($arr[0] eq "Sample"){
		print OUTPUT "$_\tUse_Gene_Name_Modify\tUse_Protein_Name_Modify\n";
	}else{
		if($arr[41] ne "RefSeq"){
			my $str=$arr[41]."\t".$arr[42]."\t".$arr[43];
			if(exists $infohash{$str}){
				print OUTPUT $_."\t".$infohash{$str}."\n";
			}else{
				print "ERROR: do not have $str in ${filename}_last_use_result_needChangeInfo.txt\n";
			}
		}else{
			if((! exists $wphash{$arr[43]})){
				my $str=$arr[41]."\t".$arr[42]."\t".$arr[43];
				if(exists $infohash{$str}){
					print OUTPUT $_."\t".$infohash{$str}."\n";
				}else{
					print "ERROR: do not have $str in ${filename}_last_use_result_needChangeInfo.txt\n";
				}
			}else{
				$arr[43]=~s/ /_/g;
				my $str2=$arr[43]."(".$arr[20].")";
				if($arr[20] eq "NA"){
					$str2=$arr[43]."(pseudo)";
				}
				print OUTPUT $_."\t".$str2."\t".$str2."\n";
			}
		}
	}
}
close INPUT;
close OUTPUT;
