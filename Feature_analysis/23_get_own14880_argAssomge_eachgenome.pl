#!/usr/bin/perl -w
use strict;
use Math::Round;
my ($argfile,$MGEmergedfile,$outfile)=@ARGV;
#argfile: ~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last.txt
#MGEmergedfile:~/MGE/MGE_result/Salmonella/he20S15/he20S15_MGE_merged_site.txt
#outfile: ~/MGE/MGE_result/Salmonella/he20S15/he20S15_MGE_merged_ARGcombine_result.txt

my @arr0=split(/\//,$MGEmergedfile);
my $sample=$arr0[$#arr0]=~s/_MGE_merged_site.txt//r;


#my %argfile_hash;
my %card_hash;
my %card_hash_extended;
my %card_hash_extended_pre;
open INPUT1,$argfile or die "Cannot open $argfile";
while(<INPUT1>){
	next if(/^ORF_ID/);
        chomp;
        my @arr=split(/\t/);
	next if($arr[0] eq "Sample");
	next if($arr[0] ne $sample);
	my $start=$arr[4];
	my $end=$arr[5];
	my $strand=$arr[6];
	my $symbol=$arr[2];
	my $id=$arr[1];
	if($arr[4] > $arr[5]){
		print "ERROR in $argfile, start not before end.\n";
                $start=$arr[5];
                $end=$arr[4];
        }
	my $str22=$symbol."\t".$start."\t".$end."\t".$strand;
	#$argfile_hash{$str22}=$arr[10]."\t".$arr[0]."\t".$arr[6]."\t".$arr[7]."\t".$arr[8]."\t".$arr[9]."\t".$arr[11]."\t".$arr[14]."\t".$arr[15]."\t".$arr[16];
	if(exists $card_hash{$symbol}{$start}){
		$card_hash{$symbol}{$start}=$card_hash{$symbol}{$start}.";;".$start."\t".$end."\t".$strand."\t".$id;
        }else{
		$card_hash{$symbol}{$start}=$start."\t".$end."\t".$strand."\t".$id;
        }
	my $start_extended=$start-10000;
	if($start_extended<0){
		$start_extended=1;
	}
	my $end_extended=$end+10000;
	if(exists $card_hash_extended{$symbol}{$start_extended}){
		$card_hash_extended{$symbol}{$start_extended}=$card_hash_extended{$symbol}{$start_extended}.";;".$start_extended."\t".$end_extended."\t".$strand."\t".$id;
		$card_hash_extended_pre{$symbol}{$start_extended}=$card_hash_extended_pre{$symbol}{$start_extended}.";;".$start."\t".$end."\t".$strand."\t".$id;
	}else{
		$card_hash_extended{$symbol}{$start_extended}=$start_extended."\t".$end_extended."\t".$strand."\t".$id;
		$card_hash_extended_pre{$symbol}{$start_extended}=$start."\t".$end."\t".$strand."\t".$id;
	}
}
close INPUT1;

## All mge do not use strand compare to cds;
my %is1;
my %is2;
my %is3;
# ISESccan results are in %is1, DANMEL results are in %is2, MobileElementFinder results are in %is3.
open INPUT2_0,$MGEmergedfile or die "Cannot open $MGEmergedfile\n";
while(<INPUT2_0>){
	chomp;
	next if(/^Sample/);
	my @arr=split(/\t/);
	next if($arr[2] ne "IS");
	if($arr[1] eq "ISEScan"){
		next if($arr[9] eq "new");
		if(! exists $is1{$arr[3]}{$arr[9]}){
			$is1{$arr[3]}{$arr[9]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[9];
		}else{
			$is1{$arr[3]}{$arr[9]}=$is1{$arr[3]}{$arr[9]}.";;;".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[9];
		}
	}elsif($arr[1] eq "DANMEL"){
			if($arr[10]!~/;/){
				if(! exists $is2{$arr[3]}{$arr[9]}){
					$is2{$arr[3]}{$arr[9]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[10];
				}else{
					$is2{$arr[3]}{$arr[9]}=$is2{$arr[3]}{$arr[9]}.";;;".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[10];
				}
			}else{
				my @tmparr=split(/;/);
				if(! exists $is2{$arr[3]}{$arr[9]}){
                                        $is2{$arr[3]}{$arr[9]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$tmparr[0];
                                }else{
                                        $is2{$arr[3]}{$arr[9]}=$is2{$arr[3]}{$arr[9]}.";;;".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$tmparr[0];
                                }
			}
	}else{
		if($arr[9] ne "."){
			if($arr[10]!~/;/){
				if(! exists $is3{$arr[3]}{$arr[9]}){
		                        $is3{$arr[3]}{$arr[9]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[10];
        		        }else{  
                		        $is3{$arr[3]}{$arr[9]}=$is3{$arr[3]}{$arr[9]}.";;;".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[10];
                		}
			}else{
				my @tmparr=split(/;/);
                                if(! exists $is3{$arr[3]}{$arr[9]}){
                                        $is3{$arr[3]}{$arr[9]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$tmparr[0];
                                }else{
                                        $is3{$arr[3]}{$arr[9]}=$is3{$arr[3]}{$arr[9]}.";;;".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$tmparr[0];
                                }
			}
		}else{
			if($arr[10]=~/;/){
				my @tmparr=split(/;/);
	                        if(! exists $is3{$arr[3]}{$tmparr[0]}){
        	                        $is3{$arr[3]}{$tmparr[0]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$tmparr[0];
                	        }else{  
                        	        $is3{$arr[3]}{$tmparr[0]}=$is3{$arr[3]}{$tmparr[0]}.";;;".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$tmparr[0];
                        	}
                        }else{
                                if(! exists $is3{$arr[3]}{$arr[10]}){
                                        $is3{$arr[3]}{$arr[10]}=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[10];
                                }else{
                                        $is3{$arr[3]}{$arr[10]}=$is3{$arr[3]}{$arr[10]}.";;;".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[10];
                                }
			}
		}
	}
}
close INPUT2_0;

my %card_hash_use;
my %mgehash;
my $abc=0;
my %printhash;
my $printstr;
my %is;
open INPUT2,$MGEmergedfile or die "Cannot open $MGEmergedfile\n";
open OUTPUT2,">$outfile" or die "Cannot write $outfile\n";
print OUTPUT2 "Sample\tSeqid\tMGE_Software\tMGE_Type\tMGE_Start\tMGE_End\tMGE_Strand\tMGE_Hit_Length\tMGE_Contig_Length\tMGE_Element_type\tMGE_Anno\tMGE_Accession\tMGE_Identity\tMGE_SCoverage\tMGE_QCoverage\tMGE_MergedInfo\tMGE_Transferability\tGeneid_Start\tGeneid_End\tGeneid_Strand\tGeneid_ID\tDistance\tPairedISname\n";
while(<INPUT2>){
	chomp;
	next if(/^Sample/);
	my @arr=split(/\t/);
	if(exists $card_hash{$arr[3]}){
#		if($arr[2] eq "transposon" || $arr[2] eq "plasmid" || $arr[2] eq "phage" || $arr[2] eq "ICE" || $arr[2] eq "IME"){
		if($arr[2]=~/transposon/ || $arr[2] eq "plasmid" || $arr[2] eq "phage" || $arr[2] eq "ICE" || $arr[2] eq "T4SS_ICE" || $arr[2] eq "AICE" || $arr[2] eq "integron"){
		### Find cds between mge.
			$abc=1;
			%card_hash_use=%card_hash;
			foreach my $key (sort keys %{$card_hash_use{$arr[3]}}){
				my @arr1=split(";;",$card_hash_use{$arr[3]}{$key});
				foreach my $arr1_key (@arr1){
					my @arr2=split("\t",$arr1_key);
					my $mge_start=$arr[4];
					my $mge_end=$arr[5];
					if($arr[1] eq "PlasmidFinder"){
						$mge_start=1;
						$mge_end=$arr[8];
					}
					if($mge_start>$mge_end){
                                                $mge_start=$mge_end;
                                                $mge_end=$arr[4];
                                        }
					if(($mge_start<=$arr2[0]) && ($mge_end>=$arr2[1])){
						$printstr="$arr[0]\t$arr[3]\t$arr[1]\t$arr[2]\t$mge_start\t$mge_end\t$arr[6]\t$arr[7]\t$arr[8]\t$arr[9]\t$arr[10]\t$arr[11]\t$arr[12]\t$arr[13]\t$arr[14]\t$arr[15]\t$arr[16]\t$arr2[0]\t$arr2[1]\t$arr2[2]\t$arr2[3]\t0\t.\n";
						$printhash{$printstr}=1;
                                                print OUTPUT2 $printstr;
					}else{next;}
				}
			}
		}
		if($arr[2] eq "integron"){
		### Find mge within the range of 20 kb in the upstream and downstream of cds area (upstream 10kb downstream 10kb). No need for mge to be completely in the region, only if have one base in.
			$abc=1;
			%card_hash_use=%card_hash_extended;
			foreach my $key (sort keys %{$card_hash_use{$arr[3]}}){
                                my @arr1=split(";;",$card_hash_use{$arr[3]}{$key});
				my @arr1_pre=split(";;",$card_hash_extended_pre{$arr[3]}{$key});
				my $tempi=0;
                                foreach my $arr1_key (@arr1){
					$tempi++;
                                        my @arr2=split("\t",$arr1_key);
					my @arr2_pre=split("\t",$arr1_pre[$tempi-1]);
					my $mge_start=$arr[4];
                                        my $mge_end=$arr[5];
                                        if($mge_start>$mge_end){
                                                $mge_start=$mge_end;
                                                $mge_end=$arr[4];
                                        }
					if($mge_start>=$arr2_pre[1] && $mge_start>=$arr2[0] && $mge_start<=$arr2[1]){
						my $distance=$mge_start-$arr2_pre[1];
						$printstr="$arr[0]\t$arr[3]\t$arr[1]\t$arr[2]\t$mge_start\t$mge_end\t$arr[6]\t$arr[7]\t$arr[8]\t$arr[9]\t$arr[10]\t$arr[11]\t$arr[12]\t$arr[13]\t$arr[14]\t$arr[15]\t$arr[16]\t$arr2_pre[0]\t$arr2_pre[1]\t$arr2_pre[2]\t$arr2_pre[3]\t$distance\t.\n";
						if(! exists $printhash{$printstr}){
							print OUTPUT2 $printstr;
						}
					}elsif($mge_end<=$arr2_pre[0] && $mge_end>=$arr2[0] && $mge_end<=$arr2[1]){
						my $distance=$arr2_pre[0]-$mge_end;
                                                $printstr="$arr[0]\t$arr[3]\t$arr[1]\t$arr[2]\t$mge_start\t$mge_end\t$arr[6]\t$arr[7]\t$arr[8]\t$arr[9]\t$arr[10]\t$arr[11]\t$arr[12]\t$arr[13]\t$arr[14]\t$arr[15]\t$arr[16]\t$arr2_pre[0]\t$arr2_pre[1]\t$arr2_pre[2]\t$arr2_pre[3]\t$distance\t.\n";
                                                if(! exists $printhash{$printstr}){
                                                        print OUTPUT2 $printstr;
                                                }
					}else{next;}
				}
			}
		}
		if($arr[2] eq "IS"){
			next if($arr[1] eq "ISEScan" && $arr[9] eq "new");
                        $abc=1;
			my $name=$arr[9];
			my $name10="";
                        if($arr[10]=~/;/){
                                my @tmparrname=split(/;/);
                                $name10=$tmparrname[0];
                        }else{
                                $name10=$arr[10];
                        }

			if($arr[1] eq "ISEScan"){
				%is=%is1;
				$name10=$arr[9];
			}elsif($arr[1] eq "DANMEL"){
				%is=%is2;
			}else{
				%is=%is3;
				if($arr[9] eq "."){
					if($arr[10]=~/;/){
	                                        my @tmparrname=split(/;/);
        	                                $name=$tmparrname[0];
                	                }else{
                        	                $name=$arr[10];
                                	}
				}
			}
			if(exists $is{$arr[3]}{$name} && $is{$arr[3]}{$name}=~/;;;/){
		                %card_hash_use=%card_hash_extended;
	        	        foreach my $key (sort keys %{$card_hash_use{$arr[3]}}){
                       		        my @arr1=split(";;",$card_hash_use{$arr[3]}{$key});
               	        	        my @arr1_pre=split(";;",$card_hash_extended_pre{$arr[3]}{$key});
       	                        	my $tempi=0;
	                                foreach my $arr1_key (@arr1){
       	                        	        $tempi++;
               	        	                my @arr2=split("\t",$arr1_key);
               		                        my @arr2_pre=split("\t",$arr1_pre[$tempi-1]);
       	                	                my $mge_start=$arr[4];
                                	        my $mge_end=$arr[5];
                                       		if($mge_start>$mge_end){
                               	                	$mge_start=$mge_end;
                        	                        $mge_end=$arr[4];
       	        	                        }
       		                                if($mge_start>=$arr2_pre[1] && $mge_start>=$arr2[0] && $mge_start<=$arr2[1]){
                	                                my @temparr0=split(/;;;/,$is{$arr[3]}{$name});
							my $mgestr=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$name10;
							my $distance1=10000;
							my $pairedis=0;
							my $pairedisname=".";
							foreach my $temparr0_0(@temparr0){
								next if($temparr0_0 eq $mgestr);
								my @temparr1=split(/\t/,$temparr0_0);
								if($temparr1[1]<=$arr2_pre[0] && $temparr1[1]>=$arr2[0] && $temparr1[1]<=$arr2[1]){
									$pairedis=1;
									my $distance1_0=$arr2_pre[0]-$temparr1[1];
									if($distance1_0 <= $distance1){
										$distance1=$distance1_0;
										$pairedisname=$temparr1[3];
									}
								}
							}
							if($pairedis>0){
								##print distance as distance1,distance2. distance1 means mge at the left of aro, distance2 means mge at the right of aro.
								my $distance=$mge_start-$arr2_pre[1];
								$printstr="$arr[0]\t$arr[3]\t$arr[1]\t$arr[2]\t$mge_start\t$mge_end\t$arr[6]\t$arr[7]\t$arr[8]\t$arr[9]\t$arr[10]\t$arr[11]\t$arr[12]\t$arr[13]\t$arr[14]\t$arr[15]\t$arr[16]\t$arr2_pre[0]\t$arr2_pre[1]\t$arr2_pre[2]\t$arr2_pre[3]\t${distance1},${distance}\t$pairedisname\n";
								if(! exists $printhash{$printstr}){
									print OUTPUT2 $printstr;
								}
							}
        	                                }elsif($mge_end<=$arr2_pre[0] && $mge_end>=$arr2[0] && $mge_end<=$arr2[1]){
							my @temparr0=split(/;;;/,$is{$arr[3]}{$name});
               	                                        my $mgestr=$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$name10;
                       	                                my $distance1=10000;
                               	                        my $pairedis=0;
							my $pairedisname=".";
                                       	                foreach my $temparr0_0(@temparr0){
                                               	                next if($temparr0_0 eq $mgestr);
                                                       	        my @temparr1=split(/\t/,$temparr0_0);
                                                               	if($temparr1[0]>=$arr2_pre[1] && $temparr1[0]>=$arr2[0] && $temparr1[0]<=$arr2[1]){
                                                                       	$pairedis=1;
                                                                        my $distance1_0=$temparr1[0]-$arr2_pre[1];
       	                                                                if($distance1_0 <= $distance1){
               	                                                                $distance1=$distance1_0;
										$pairedisname=$temparr1[3];
                       	                                                }
                               	                                }
                                       	                }
                                               	        if($pairedis>0){
	                                        	        my $distance=$arr2_pre[0]-$mge_end;
       	                                        		$printstr="$arr[0]\t$arr[3]\t$arr[1]\t$arr[2]\t$mge_start\t$mge_end\t$arr[6]\t$arr[7]\t$arr[8]\t$arr[9]\t$arr[10]\t$arr[11]\t$arr[12]\t$arr[13]\t$arr[14]\t$arr[15]\t$arr[16]\t$arr2_pre[0]\t$arr2_pre[1]\t$arr2_pre[2]\t$arr2_pre[3]\t${distance},${distance1}\t$pairedisname\n";
                	                        	        if(! exists $printhash{$printstr}){
       	                	        	                        print OUTPUT2 $printstr;
               	        		                        }
               		                	        }else{next;}
       	                	        	}
					}
         	               }
			}else{
				next;
			}
		}
		if($abc==0){
			print "ERROR: $sample in merge cds, $arr[2] in MGE_merged_site.txt is not element we want.\n";
		}
	}
	
}
close INPUT2;
close OUTPUT2;
