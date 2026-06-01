use strict;
use warnings;

my $outfile=$ARGV[0];

my %refseqresult;
open INPUT0,"/lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno/RefSeq_Anno_Result_14888file_addinfo.txt" or die "cannot open RefSeq_Anno_Result_14888file_addinfo.txt\n";
while(<INPUT0>){
	chomp;
	my @arr=split(/\t/);
	next if($arr[0] eq "qseqid");
	my @arr2=split(/ /,$arr[13]);
	my $rsgene="";
	my $rslocustag="";
	my $rsprotein="";
	my $rsproteinid="";
	my $rsgbkey="";
	if($arr[13]=~/\[gene=(.*?)\]/){
		$rsgene=$1;
	}else{
		$rsgene="NA";
	}
	if($arr[13]=~/\[locus_tag=(.*?)\]/){
                $rslocustag=$1;
        }else{
                $rslocustag="NA";
        }
	if($arr[13]=~/\[protein=(.*?)\]/){
                $rsprotein=$1;
        }else{
                $rsprotein="NA";
        }
	if($arr[13]=~/\[protein_id=(.*?)\]/){
                $rsproteinid=$1;
        }else{
                $rsproteinid="NA";
        }
	if($arr[13]=~/\[gbkey=(.*?)\]/){
                $rsgbkey=$1;
        }else{
                $rsgbkey="NA";
        }
	
	$refseqresult{$arr[0]}{$arr[1]}=$arr[5]."\t".$arr[6]."\t".$arr[7]."\t".$arr[8]."\t".$arr[9]."\t".$arr[10]."\t".$arr[11]."\t".$arr[12]."\t".$arr[13]."\t".$rsgene."\t".$rslocustag."\t".$rsprotein."\t".$rsproteinid."\t".$rsgbkey;
}
close INPUT0;


open OUTPUT,">$outfile" or die $!;
print OUTPUT "Sample\tGene_id\tGene_contig\tGene_type\tGene_start\tGene_end\tGene_strand\tGene_phase\tRefSeq_alignment_length\tRefSeq_mismatches\tRefSeq_gap_opens\tRefSeq_identity\tRefSeq_qcov\tRefSeq_scov\tRefSeq_evalue\tRefSeq_bit_score\tRefSeq_subject_id\tRefSeq_gene\tRefSeq_locus_tag\tRefSeq_protein\tRefSeq_protein_id\tRefSeq_gbkey\tProkka_gene\tProkka_product\tARGblastn_alignment_length\tARGblastn_mismatches\tARGblastn_gap_opens\tARGblastn_identity\tARGblastn_qcov\tARGblastn_scov\tARGblastn_evalue\tARGblastn_bit_score\tARGblastn_subject_id\tARGblastn_gene_name\tARGblastn_alternative_name\tARGblastn_product_name\tARGblastn_gene_family\tARGblastn_drug_class\tARGblastn_drug_name\tARGblastn_mechanism_of_resistance\tARGblastn_other_notes\n";
open INPUT,"/lustre/home/niejingyi2023/project/bjCDC/Salmonella_MIC_usegenome.list" or die $!;
while(my $sample=<INPUT>){
	chomp($sample);
	my %genes=();
	my %prokkaresult=();
	open INPUT1,"/lustre/home/niejingyi2023/project/bjCDC/MGE/MGE_result/Salmonella/$sample/prokka/$sample.gff" or die "cannot open MGE/MGE_result/Salmonella/$sample/prokka/$sample.gff\n";
	while(<INPUT1>){
		chomp;
		next if(/^#/);
		if(/^>/){last;}
		my @arr=split(/\t/);
		my $prokkageneid="";
		my $prokkagene="";
                my $prokkaprotein="";
		if($arr[8]=~/ID=(.*?);/){
			$prokkageneid=$1;
			$genes{$prokkageneid}=$sample."\t".$prokkageneid."\t".$arr[0]."\t".$arr[2]."\t".$arr[3]."\t".$arr[4]."\t".$arr[6]."\t".$arr[7];
		}
		if($arr[8]=~/gene=(.*?);/){
			$prokkagene=$1;
                }else{
			$prokkagene="NA";
		}
		if($arr[8]=~/product=(.*)/){
                        $prokkaprotein=$1;
                }else{
                        $prokkaprotein="NA";
                }

		$prokkaresult{$prokkageneid}=$prokkagene."\t".$prokkaprotein;
		
	}
	close INPUT1;

	my %argdiamondresult=();
	my $tail="diamond.out.processed.txt";
	if($outfile=~/add_mutation_models/){
		my $tail="diamond.out_addpoint_models.processed.txt";	
	}
	open INPUT2,"/lustre/home/niejingyi2023/project/bjCDC/MGE/MGE_result/Salmonella/$sample/ARG_combine_diamond/$sample.$tail" or die "cannot open MGE/MGE_result/Salmonella/$sample/ARG_combine_diamond/$sample.$tail\n";
	while(<INPUT2>){
		chomp;
		my @arr=split(/\t/);
		next if($arr[0] eq "query_id");
		$argdiamondresult{$arr[0]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[11]."\t".$arr[12]."\t".$arr[13]."\t".$arr[14]."\t".$arr[15]."\t".$arr[16]."\t".$arr[17]."\t".$arr[18]."\t".$arr[19]."\t".$arr[20]."\t".$arr[21]."\t".$arr[22]."\t".$arr[23]."\t".$arr[24];
	}
	close INPUT2;


	foreach my $gene(sort keys %genes){
		print OUTPUT $genes{$gene};
		if(exists $refseqresult{$sample}{$gene}){
			print OUTPUT "\t".$refseqresult{$sample}{$gene};
		}else{
			print OUTPUT "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		}
		if(exists $prokkaresult{$gene}){
			print OUTPUT "\t".$prokkaresult{$gene};
		}else{
			print "ERROR: do not have $gene in prokka file\n";
		}
		if(exists $argdiamondresult{$gene}){
                        print OUTPUT "\t".$argdiamondresult{$gene}."\n";
                }else{
                        print OUTPUT "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
                }
	}
}
close INPUT;
close OUTPUT;
