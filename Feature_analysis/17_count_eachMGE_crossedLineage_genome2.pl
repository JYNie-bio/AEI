use strict;
use warnings;

my ($input0,$input,$output)=(@ARGV);

my %hash;
open INPUT0,"$input0" or die $!;
while(<INPUT0>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "Lineage_level");
	$hash{$arr[0]}=$arr[1];
}
close INPUT0;

#my %hash_kingdom;
my %hash_phylum;
my %hash_class;
my %hash_order;
my %hash_family;
my %hash_genus;
my %mgetype;
open INPUT1,"$input" or die $!;
while(<INPUT1>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "MGE_Type");
        my $str="$arr[0]\t$arr[1]";
	if(! exists $mgetype{$str}){
		$mgetype{$str}=$arr[2];
	}elsif($mgetype{$str} !~ /\b\Q$arr[2]\E\b/){
		$mgetype{$str}=$mgetype{$str}.";".$arr[2];
	}
#	if(! exists $hash_kingdom{$str}){
#		$hash_kingdom{$str}=$arr[7];
#	}elsif($hash_kingdom{$str} !~ /\b\Q$arr[7]\E\b/){
#		$hash_kingdom{$str}=$hash_kingdom{$str}.";".$arr[7];
#	}
	if(! exists $hash_phylum{$str}){
		$hash_phylum{$str}=$arr[8];
        }elsif($hash_phylum{$str} !~ /\b\Q$arr[8]\E\b/){
		$hash_phylum{$str}=$hash_phylum{$str}.";".$arr[8];
	}
	if(! exists $hash_class{$str}){
		$hash_class{$str}=$arr[9];
        }elsif($hash_class{$str} !~ /\b\Q$arr[9]\E\b/){
		$hash_class{$str}=$hash_class{$str}.";".$arr[9];
	}
	if(! exists $hash_order{$str}){
		$hash_order{$str}=$arr[10];
        }elsif($hash_order{$str} !~ /\b\Q$arr[10]\E\b/){
		$hash_order{$str}=$hash_order{$str}.";".$arr[10];
	}
	if(! exists $hash_family{$str}){
		$hash_family{$str}=$arr[11];
	}elsif($hash_family{$str} !~ /\b\Q$arr[11]\E\b/){
		$hash_family{$str}=$hash_family{$str}.";".$arr[11];
	}
	if(! exists $hash_genus{$str}){
		$hash_genus{$str}=$arr[12];
        }elsif($hash_genus{$str} !~ /\b\Q$arr[12]\E\b/){
		$hash_genus{$str}=$hash_genus{$str}.";".$arr[12];
	}
}
close INPUT1;

open OUTPUT,">$output" or die $!;
print OUTPUT "MGE_Type\tElement_Type\tCrossedLineage\tCrossedName\tCrossedLineageNumber\tLineageAllNumber\tCrossedLineageProportion\n";
foreach my $mge(sort keys %mgetype){
#	my @kingdomarr=split(/;/,$hash_kingdom{$mge});
#	my $kingdomnum=$#kingdomarr+1;
#	if($kingdomnum>=2){
#		print OUTPUT $mge."\tKingdom\t".$hash_kingdom{$mge}."\n";
#	}

                	my @phylumarr=split(/;/,$hash_phylum{$mge});
                	my $phylumnum=$#phylumarr+1;
			my $phylumprop=$phylumnum/$hash{"Phylum"};
                	if($phylumnum>=2 || ($phylumnum==1 && ($hash_phylum{$mge} ne "Proteobacteria"))){
                       		print OUTPUT $mge."\tPhylum\t".$hash_phylum{$mge}."\t".$phylumnum."\t".$hash{"Phylum"}."\t".$phylumprop."\n";
			}
                        my @classarr=split(/;/,$hash_class{$mge});
                        my $classnum=$#classarr+1;
			my $classprop=$classnum/$hash{"Class"};
                        if($classnum>=2 || ($classnum==1 && ($hash_class{$mge} ne "Gammaproteobacteria"))){
                                print OUTPUT $mge."\tClass\t".$hash_class{$mge}."\t".$classnum."\t".$hash{"Class"}."\t".$classprop."\n";
                        }
                        my @orderarr=split(/;/,$hash_order{$mge});
                        my $ordernum=$#orderarr+1;
			my $orderprop=$ordernum/$hash{"Order"};
                        if($ordernum>=2 || ($ordernum==1 && ($hash_order{$mge} ne "Enterobacterales"))){
                                print OUTPUT $mge."\tOrder\t".$hash_order{$mge}."\t".$ordernum."\t".$hash{"Order"}."\t".$orderprop."\n";
                        }
                        my @familyarr=split(/;/,$hash_family{$mge});
                        my $familynum=$#familyarr+1;
			my $familyprop=$familynum/$hash{"Family"};
                        if($familynum>=2 || ($familynum==1 && ($hash_family{$mge} ne "Enterobacteriaceae"))){
                                print OUTPUT $mge."\tFamily\t".$hash_family{$mge}."\t".$familynum."\t".$hash{"Family"}."\t".$familyprop."\n";
                        }
                        my @genusarr=split(/;/,$hash_genus{$mge});
                        my $genusnum=$#genusarr+1;
			my $genusprop=$genusnum/$hash{"Genus"};
                        if($genusnum>=2 || ($genusnum==1 && ($hash_genus{$mge} ne "Salmonella"))){
                                print OUTPUT $mge."\tGenus\t".$hash_genus{$mge}."\t".$genusnum."\t".$hash{"Genus"}."\t".$genusprop."\n";
                        }
			my @speciesarr=split(/;/,$mgetype{$mge});
			my $speciesnum=$#speciesarr+1;
			my $speciesprop=$speciesnum/$hash{"Species"};
			if($speciesnum>=2 || ($speciesnum==1 && ($mgetype{$mge} ne "Salmonella_enterica" || $mgetype{$mge} ne "Salmonella_sp."))){
				print OUTPUT $mge."\tSpecies\t".$mgetype{$mge}."\t".$speciesnum."\t".$hash{"Species"}."\t".$speciesprop."\n";
			}
}
close OUTPUT;
