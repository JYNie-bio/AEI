open INPUT1,"~/Feature_analysis/updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_allMGE_crossedLineage_genome1.txt" or die $!;
my %hash;
my %hashmge;
while(<INPUT1>){
        chomp;
        my @arr=split(/\t/);
        next if($arr[0] eq "MGE_Type");
	my $str=$arr[0]."\t".$arr[1];
	$hashmge{$str}=1;
	if(! exists $hash{$str}{$arr[2]}){
		$hash{$str}{$arr[2]}="$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]";
	}else{
		$hash{$str}{$arr[2]}=$hash{$str}{$arr[2]}.";;;"."$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]";
	}
}
close INPUT1;

open OUTPUT,">~/Feature_analysis/updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_allMGE_crossedLineage_genome1_lineageGenome2.txt";
print OUTPUT "MGE_Type\tElement_Type\tCrossedLineage\tCrossedName\n";
foreach my $key(sort keys %hashmge){
	my @lineages=keys %{$hash{$key}};
	if(grep /Phylum/, @lineages){
		my @arr=split(/;;;/,$hash{$key}{"Phylum"});
		foreach my $line(@arr){
			print OUTPUT $line."\n";
		}
	}elsif(grep /Class/, @lineages){
		my @arr=split(/;;;/,$hash{$key}{"Class"});
                foreach my $line(@arr){
                        print OUTPUT $line."\n";
                }
	}elsif(grep /Order/, @lineages){
                my @arr=split(/;;;/,$hash{$key}{"Order"});
                foreach my $line(@arr){
                        print OUTPUT $line."\n";
                }
        }elsif(grep /Family/, @lineages){
                my @arr=split(/;;;/,$hash{$key}{"Family"});
                foreach my $line(@arr){
                        print OUTPUT $line."\n";
                }
        }elsif(grep /Genus/, @lineages){
                my @arr=split(/;;;/,$hash{$key}{"Genus"});
                foreach my $line(@arr){
                        print OUTPUT $line."\n";
                }
        }else{
		my @arr=split(/;;;/,$hash{$key}{"Species"});
                foreach my $line(@arr){
                        print OUTPUT $line."\n";
                }
	}
}
close OUTPUT;
