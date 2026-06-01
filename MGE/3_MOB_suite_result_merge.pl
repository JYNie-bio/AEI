use strict;
use warnings;

open OUTPUT,">~/MGE/MGE_result/Salmonella_MOBsuite/14880_MOB_suite_result.txt";
print OUTPUT "genome\tsample_id\tnum_contigs\tsize\tgc\tmd5\trep_type(s)\trep_type_accession(s)\trelaxase_type(s)\trelaxase_type_accession(s)\tmpf_type\tmpf_type_accession(s)\torit_type(s)\torit_accession(s)\tpredicted_mobility\tmash_nearest_neighbor\tmash_neighbor_distance\tmash_neighbor_identification\tprimary_cluster_id\tsecondary_cluster_id\tpredicted_host_range_overall_rank\tpredicted_host_range_overall_name\tobserved_host_range_ncbi_rank\tobserved_host_range_ncbi_name\treported_host_range_lit_rank\treported_host_range_lit_name\tassociated_pmid(s)\n";

my $mgepath="~/MGE/MGE_result/Salmonella_MOBsuite/MOBsuite_result";
opendir DIR,"$mgepath" or die $!;
my @files=readdir(DIR);
closedir DIR;
foreach my $file(@files){
	next if($file eq ".");
        next if($file eq "..");
	my $str=$file=~s/_mobsuite.txt//gr;
		open INPUT,"$mgepath/$file";
		while(<INPUT>){
			my @arr=split(/\t/);
                        next if($arr[0] eq "sample_id");
			if($#arr>3){
				print OUTPUT $str."\t".$_;
			}
                }
                close INPUT;
}
close OUTPUT;

