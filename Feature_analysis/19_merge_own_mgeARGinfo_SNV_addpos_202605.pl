use strict;
use warnings;

my $temp="";
if($ARGV[0]=="1"){
        $temp="_add_mutation_models";
}

my %snps;
open INPUT,"~/Feature_analysis/use_topfeature_name${temp}.txt" or die $!;
while(<INPUT>){
        chomp;
	my @arr=split(/\t/);
        $snps{$arr[0]}=1;
}
close INPUT;


opendir DIR,"~/ARGSNP/mummer${temp}";
my @files0=readdir(DIR);
closedir DIR;

my %hashprotein;
my %hash;
my %samples;
foreach my $file0(@files0){
        next if($file0!~/snps.snp.indel959090.filtedARG$/);
        my $sample=$file0=~s/.snps.snp.indel959090.filtedARG//r;
        $samples{$sample}=1;
        open INPUT2,"~/ARGSNP/mummer${temp}/$file0" or die "cannot open ~/ARGSNP/mummer${temp}/$file0\n";
        while(<INPUT2>){
                chomp;
                my @arr=split(/\t/);
                next if($arr[0] eq "Ref_contig");
                my @arr2=split(/\|/,$arr[0]);
		$arr2[4]=~s/cryptic_aminoglycoside_nucleotidyltransferase_ANT\(3''\)\/ANT\(9\)/ANT(3'')-Ia_family_aminoglycoside_nucleotidyltransferase_AadA/;
                $arr2[4]=~s/aminoglycoside_6'-N-acetyltransferase_AAC\(6'\)-Iy/aminoglycoside_N-acetyltransferase_AAC(6')-Iy/;
                my $str=$arr2[4]."|".$arr[3]."-".$arr[1]."-".$arr[2];
		next if(! exists $snps{$str});
                $hash{$sample}{$str}=1;
		if(! exists $hashprotein{$sample}{$arr[9]}){
                        $hashprotein{$sample}{$arr[9]}=$str;
                }else{
                        if($hashprotein{$sample}{$arr[9]} !~ /\b\Q$str\E\b/){
                                $hashprotein{$sample}{$arr[9]}=$hashprotein{$sample}{$arr[9]}.";".$str;
                        }
                }
        }
        close INPUT2;
}


my %hashmge;
my %genomes;
my $mgepath="~/MGE/MGE_result/Salmonella";
	foreach my $file(sort keys %hashprotein){
                $genomes{$file}=1;
		my $str;
		my %hashplasmid=();
		open INPUT,"$mgepath/$file/${file}_MGE_merged_ARGcombine_result${temp}_202605.txt" or die "cannot open $mgepath/$file/${file}_MGE_merged_ARGcombine_result${temp}_202605.txt\n";
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
		open INPUT,"$mgepath/$file/${file}_MGE_merged_ARGcombine_result${temp}_202605.txt" or die "cannot open $mgepath/$file/${file}_MGE_merged_ARGcombine_result${temp}_202605.txt\n";
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
				my @temparr=split(/;/,$hashprotein{$arr[0]}{$arr[20]});
                                foreach my $tempfeature(@temparr){
                                        if(! exists $hashmge{$arr[0]}{$tempfeature}){
                                                $hashmge{$arr[0]}{$tempfeature}=$str;
                                        }else{
                                                if($hashmge{$arr[0]}{$tempfeature} !~ /\b\Q$str\E\b/){
                                                        $hashmge{$arr[0]}{$tempfeature}=$hashmge{$arr[0]}{$tempfeature}.";;".$str;
                                                }
                                        }
                                }
			}elsif($arr[2] eq "BacAnt"){
                                $str=$arr[3]."\t".$arr[9];
				if(exists $hashplasmid{$arr[1]}){
                                        $str=$str."\tPlasmid";
                                }else{
                                        $str=$str."\tChromosome";
                                }
				my @temparr=split(/;/,$hashprotein{$arr[0]}{$arr[20]});
                                foreach my $tempfeature(@temparr){
                                        if(! exists $hashmge{$arr[0]}{$tempfeature}){
                                                $hashmge{$arr[0]}{$tempfeature}=$str;
                                        }else{
                                                if($hashmge{$arr[0]}{$tempfeature} !~ /\b\Q$str\E\b/){
                                                        $hashmge{$arr[0]}{$tempfeature}=$hashmge{$arr[0]}{$tempfeature}.";;".$str;
                                                }
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
				my @temparr=split(/;/,$hashprotein{$arr[0]}{$arr[20]});
                                foreach my $tempfeature(@temparr){
                                        if(! exists $hashmge{$arr[0]}{$tempfeature}){
                                                $hashmge{$arr[0]}{$tempfeature}=$str;
                                        }else{
                                                if($hashmge{$arr[0]}{$tempfeature} !~ /\b\Q$str\E\b/){
                                                        $hashmge{$arr[0]}{$tempfeature}=$hashmge{$arr[0]}{$tempfeature}.";;".$str;
                                                }
                                        }
                                }
                        }
		}
		close INPUT;
	}


open OUTPUT,">~/Feature_analysis/Salmonella_own_14880_mgeARGinfo_SNV${temp}_addpos0_202605.txt";
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

open INPUT,"~/Feature_analysis/Salmonella_own_14880_mgeARGinfo_SNV${temp}_addpos0_202605.txt";
open OUTPUT,">~/Feature_analysis/Salmonella_own_14880_mgeARGinfo_SNV${temp}_addpos_202605.txt";
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

