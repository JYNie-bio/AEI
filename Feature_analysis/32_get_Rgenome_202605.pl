use strict;
use warnings;

my ($antibiotic)=(@ARGV);
#antibiotic: AMP\AMS\...
#metafile: predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_year_host_country.txt
##mgefile: predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_mge.txt
##outmetafile: predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_year_host_country_onlyR.txt
##outmgefile: predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_mge_onlyR.txt

my %hashgenome;
my %hash;
open INPUT,"/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/predict_use_data_add_mutation_models/Salmonella_enterica_ownuse14880_public_184252_Allantibiotic_topfeature_matrix_result.txt";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	if($arr[0] eq "Genome"){
		for(my $i=1;$i<=$#arr;$i++){
			$hashgenome{$i}=$arr[$i];
		}
	}else{
		for(my $i=1;$i<=$#arr;$i++){
			if($arr[$i] eq "R"){
				$hash{$arr[0]}{$hashgenome{$i}}=1;
			}
		}
	}
}
close INPUT;

open INPUT1,"/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/predict_use_data_add_mutation_models/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_year_host_country.txt";
open OUTPUT1,">/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/predict_use_data_add_mutation_models/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_year_host_country_onlyR.txt";
while(<INPUT1>){
	my @arr=split(/\t/);
	if($arr[0] eq "Genome"){
		print OUTPUT1 $_;
	}elsif(exists $hash{$arr[0]}{$antibiotic}){
		print OUTPUT1 $_;
	}
}
close INPUT1;
close OUTPUT1;

open INPUT1,"/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/predict_use_data_add_mutation_models/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_mge_202605.txt";
open OUTPUT1,">/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/predict_use_data_add_mutation_models/Salmonella_enterica_ownuse14846_publicANI96_173348_${antibiotic}_topfeature_matrix_add_mge_onlyR_202605.txt";
while(<INPUT1>){
        my @arr=split(/\t/);
        if($arr[0] eq "Genome"){
                print OUTPUT1 $_;
        }elsif(exists $hash{$arr[0]}{$antibiotic}){
                print OUTPUT1 $_;
        }
}
close INPUT1;
close OUTPUT1;

