my %hash;
my %hash2;
my %nohash;
my $filename=$ARGV[0];
open INPUT0,"~/RefSeq_Anno/${filename}_notPASSarg908080_PASSrefseq908080_argAssoProtein_use.txt";
while(<INPUT0>){
        chomp;
	my @arr=split(/\t/);
	if($#arr>0){
		my @arr2=split(/;/,$arr[1]);
		foreach my $key(@arr2){
			$hash2{$arr[0]}{$key}=1;
		}
		if($#arr>1){
			$nohash{$arr[0]}{$arr[2]}=1;
		}
	}else{
		$hash{$arr[0]}=1;
	}
}
close INPUT0;


open INPUT,"~/RefSeq_Anno/${filename}_notPASSarg908080_PASSrefseq908080.txt";
open OUTPUT,">~/RefSeq_Anno/${filename}_notPASSarg908080_PASSrefseq908080_use.txt";
print OUTPUT "Sample\tGene_id\tGene_contig\tGene_type\tGene_start\tGene_end\tGene_strand\tGene_phase\tRefSeq_alignment_length\tRefSeq_mismatches\tRefSeq_gap_opens\tRefSeq_identity\tRefSeq_qcov\tRefSeq_scov\tRefSeq_evalue\tRefSeq_bit_score\tRefSeq_subject_id\tRefSeq_gene\tRefSeq_locus_tag\tRefSeq_protein\tRefSeq_protein_id\tRefSeq_gbkey\tProkka_gene\tProkka_product\tARGblastn_alignment_length\tARGblastn_mismatches\tARGblastn_gap_opens\tARGblastn_identity\tARGblastn_qcov\tARGblastn_scov\tARGblastn_evalue\tARGblastn_bit_score\tARGblastn_subject_id\tARGblastn_gene_name\tARGblastn_alternative_name\tARGblastn_product_name\tARGblastn_gene_family\tARGblastn_drug_class\tARGblastn_drug_name\tARGblastn_mechanism_of_resistance\tARGblastn_other_notes\tAnnotation_source\tUse_Gene_name\tUse_Protein_name\n";
while(<INPUT>){
	chomp;
	my @arr=split(/\t/);
	if(((exists $hash{$arr[19]}) || (exists $hash2{$arr[19]}{$arr[17]})) && (! exists $nohash{$arr[19]}{$arr[17]})){
		print OUTPUT $_."\tRefSeq\t".$arr[17]."\t".$arr[19]."\n";
	}
}
close INPUT;
