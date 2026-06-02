use strict;
use warnings;

my $temp="";
if($ARGV[0]=="1"){
        $temp="_add_mutation_models";
}

my %hashprotein;
my %hash;
open INPUT3,"~/Feature_analysis/Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_arg${temp}_last_use_result_last.txt" or die $!;
while(<INPUT3>){
        chomp;
        my @arr=split(/\t/);
        next if($arr[0] eq "Sample");
	$hash{$arr[0]}{$arr[45]}=1;
        $hashprotein{$arr[0]}{$arr[1]}=$arr[45];
}
close INPUT3;

my %hashmge;
my %genomes;
my $mgepath="~/Feature_analysis/MGE_prerunned";
	opendir DIR,"$mgepath" or die $!;
	my @files=readdir(DIR);
	closedir DIR;
	foreach my $file(@files){
		#next if($file eq ".");
		#next if($file eq "..");
		next if($file !~ /MGE_merged_ARGcombine_result${temp}_202605.txt/);
		my $temp=$file=~s/_MGE_merged_ARGcombine_result${temp}_202605.txt//r;
                $genomes{$temp}=1;
		my $str;
		my %hashplasmid=();
		open INPUT,"$mgepath/$file";
                while(<INPUT>){
                        chomp;
                        my @arr=split(/\t/);
                        next if($arr[0] eq "Sample");
                        if($#arr<10){last;}
                        next if(! exists $hashprotein{$arr[0]}{$arr[20]});
			if($arr[3] eq "plasmid"){
                                $hashplasmid{$arr[1]}=1;
                        }
		}
		close INPUT;
		open INPUT,"$mgepath/$file";
		while(<INPUT>){
			chomp;
			my @arr=split(/\t/);
			next if($arr[0] eq "Sample");
			if($#arr<10){last;}
			next if(! exists $hashprotein{$arr[0]}{$arr[20]});
			if($arr[2] eq "MobileElementFinder"){
				if($arr[3]=~/transposon/){
					$str="transposon\t".$arr[9];
				}elsif($arr[3] eq "IS"){
					$str=$arr[3]."\t".$arr[10];
				}else{
					$str=$arr[3]."\t".$arr[9];
				}
				if(exists $hashplasmid{$arr[1]}){
                                        $str=$str."\tPlasmid";
                                }else{
                                        $str=$str."\tChromosome";
                                }
				if(! exists $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}){
					$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$str;
				}else{
					if($hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}} !~ /\b\Q$str\E\b/){
						$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}.";;".$str;
					}
				}
			}elsif($arr[2] eq "BacAnt"){
                                $str=$arr[3]."\t".$arr[9];
				if(exists $hashplasmid{$arr[1]}){
                                        $str=$str."\tPlasmid";
                                }else{
                                        $str=$str."\tChromosome";
                                }
				if(! exists $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}){
                                        $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$str;
                                }else{
                                        if($hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}} !~ /\b\Q$str\E\b/){
                                                $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}.";;".$str;
                                        }
                                }
                        }elsif($arr[2] eq "plasmid"){
				if($arr[9]=~/\((.*)\)/){
					$str=$arr[3]."\t".$1;
				}else{
					$str=$arr[3]."\t".$arr[9];
				}
				if(exists $hashplasmid{$arr[1]}){
                                        $str=$str."\tPlasmid";
                                }else{
                                        $str=$str."\tChromosome";
                                }
				if(! exists $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}){
                                        $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$str;
                                }else{
                                        if($hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}} !~ /\b\Q$str\E\b/){
                                                $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}.";;".$str;
                                        }
                                }
			}
		}
		close INPUT;
	}
$mgepath="~/Feature_analysis/MGE_bc_result/Salmonella_enterica";
	opendir DIR,"$mgepath" or die $!;
	@files=readdir(DIR);
	closedir DIR;
	foreach my $file(@files){
                next if($file eq ".");
                next if($file eq "..");
		my $temp=$file=~s/_MGE_merged_ARGcombine_result${temp}_202605.txt//r;
                $genomes{$temp}=1;
		my $str;
		my %hashplasmid=();
		open INPUT,"$mgepath/$file/${file}_MGE_merged_ARGcombine_result${temp}_202605.txt";
                while(<INPUT>){
                        chomp;
                        my @arr=split(/\t/);
                        next if($arr[0] eq "Sample");
                        if($#arr<10){last;}
                        next if(! exists $hashprotein{$arr[0]}{$arr[20]});
			if($arr[3] eq "plasmid"){
                                $hashplasmid{$arr[1]}=1;
                        }
		}
		close INPUT;
		open INPUT,"$mgepath/$file/${file}_MGE_merged_ARGcombine_result${temp}_202605.txt";
		while(<INPUT>){
                        chomp;
			my @arr=split(/\t/);
                        next if($arr[0] eq "Sample");
			if($#arr<10){last;}
			next if(! exists $hashprotein{$arr[0]}{$arr[20]});
			if($arr[2] eq "MobileElementFinder"){
                                if($arr[3]=~/transposon/){
                                        $str="transposon\t".$arr[9];
                                }elsif($arr[3] eq "IS"){
                                        $str=$arr[3]."\t".$arr[10];
                                }else{  
                                        $str=$arr[3]."\t".$arr[9];
                                }
				if(exists $hashplasmid{$arr[1]}){
                                        $str=$str."\tPlasmid";
                                }else{
                                        $str=$str."\tChromosome";
                                }
                                if(! exists $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}){
                                        $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$str;
                                }else{  
                                        if($hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}} !~ /\b\Q$str\E\b/){
                                                $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}.";;".$str;
                                        }       
                                }       
                        }elsif($arr[2] eq "BacAnt"){
                                $str=$arr[3]."\t".$arr[9];
				if(exists $hashplasmid{$arr[1]}){
                                        $str=$str."\tPlasmid";
                                }else{
                                        $str=$str."\tChromosome";
                                }
                                if(! exists $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}){
                                        $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$str;
                                }else{  
                                        if($hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}} !~ /\b\Q$str\E\b/){
                                                $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}.";;".$str;
                                        }       
                                }       
                        }elsif($arr[3] eq "plasmid"){
                                if($arr[9]=~/\((.*)\)/){
                                        $str=$arr[3]."\t".$1;
                                }else{
					$str=$arr[3]."\t".$arr[9];
                                }
				if(exists $hashplasmid{$arr[1]}){
                                        $str=$str."\tPlasmid";
                                }else{
                                        $str=$str."\tChromosome";
                                }
                                if(! exists $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}){
                                        $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$str;
                                }else{
                                        if($hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}} !~ /\b\Q$str\E\b/){
                                                $hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}=$hashmge{$arr[0]}{$hashprotein{$arr[0]}{$arr[20]}}.";;".$str;
                                        }
                                }
                        }

                }
                close INPUT;
	}


open OUTPUT,">~/Feature_analysis/Salmonella_public_184251_mgeARGinfo${temp}_addpos0_202605.txt";
print OUTPUT "Genome\tFeature\tIsWith_MGE\tMGE_Type\tElement_Type\tGenome_position\n";
foreach my $sample(sort keys %hash){
	foreach my $feature(sort keys %{$hash{$sample}}){
		if(! exists $hashmge{$sample}{$feature}){
			print OUTPUT $sample."\t".$feature."\tNo\t-\t-\t-\n";
		}else{
			my @arr2=split(/;;/,$hashmge{$sample}{$feature});
				foreach my $key(@arr2){
					print OUTPUT $sample."\t".$feature."\tYes\t$key\n";
				}
		}
	}
}
close OUTPUT;

open INPUT,"~/Feature_analysis/Salmonella_public_184251_mgeARGinfo${temp}_addpos0_202605.txt";
open OUTPUT,">~/Feature_analysis/Salmonella_public_184251_mgeARGinfo${temp}_addpos_202605.txt";
my %temphash;
while(<INPUT>){
        chomp;
        if(/^Genome/){
                print OUTPUT $_."\n";
        }else{
                my @arr=split(/\t/);
                my $str="$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]";
                if(! exists $temphash{$str}){
                        $temphash{$str}=$arr[5];
                }elsif($temphash{$str} !~ /\Q$arr[5]\E/){
                        $temphash{$str}=$temphash{$str}."_".$arr[5];
                }
        }
}
close INPUT;
foreach my $key (sort keys %temphash){
        print OUTPUT $key."\t".$temphash{$key}."\n";
}
close OUTPUT;

