use strict;
use warnings;

my ($inputfile,$mgefile,$outfile)=(@ARGV);

my %hashlineage;
open INPUT1,"~/Feature_analysis/updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_allMGE_crossedLineage_genome1_lineageGenome2.txt" or die $!;
while(<INPUT1>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "MGE_Type");
	my $str=$arr[0]."\t".$arr[1];
	$hashlineage{$str}=$arr[2]."\t".$arr[3];
}
close INPUT1;

my %hashgenome0;
my %hashgenome;
my %allgenome;
open INPUT2,$inputfile or die $!;
while(<INPUT2>){
        chomp;
        my @arr=split(/\t/);
        if($arr[0] eq "Genome"){
		for(my $i=1;$i<($#arr-2);$i++){
			$hashgenome0{$i}=$arr[$i];
		}
	}else{
		for(my $i=1;$i<($#arr-2);$i++){
			$allgenome{$arr[0]}="$arr[$#arr-2]\t$arr[$#arr-1]\t$arr[$#arr]";
			if($arr[$i]>0){
				$hashgenome{$arr[0]}{$hashgenome0{$i}}="$arr[$#arr-2]\t$arr[$#arr-1]\t$arr[$#arr]";
			}
		}
	}
}
close INPUT2;


my %hashmge;
my %usegenome;
open INPUT3,$mgefile or die "$!";
open OUTPUT,">$outfile" or die $!;
print OUTPUT "Genome\tYear\tHost\tCountry\tFeature\tIsWith_MGE\tMGE_Type\tElement_Type\tCrossedLineage\tCrossedName\tFeatureMGE_genome_position\n";
while(<INPUT3>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Genome");
	if (exists $hashgenome{$arr[0]}{$arr[1]}){
		$usegenome{$arr[0]}=1;
		print OUTPUT $arr[0]."\t".$hashgenome{$arr[0]}{$arr[1]}."\t".$arr[1]."\t$arr[2]\t$arr[3]\t$arr[4]";
		my $tempstr=$arr[3]."\t".$arr[4];
		if($arr[2] eq "Yes"){
			if(exists $hashlineage{$tempstr}){
				print OUTPUT "\t".$hashlineage{$tempstr}."\t";
			}else{
				#print "ERROR: do not have $tempstr in %hashlineage.\n";
				print OUTPUT "\t-\t-\t";
			}
			print OUTPUT "$arr[5]\n";
		}else{
			print OUTPUT "\t-\t-\t-\n";
		}
	}
}
close INPUT3;
close OUTPUT;

open OUTPUT2,">>$outfile" or die $!;
foreach my $genome(sort keys %allgenome){
	if(! exists $usegenome{$genome}){
		print OUTPUT2 $genome."\t".$allgenome{$genome}."\t-\t-\t-\t-\t-\t-\t-\n";
	}
}
close OUTPUT2;
