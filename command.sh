############################### process & annotate ###################################################
cd ~/ST_Serotype
ls ~/data/Salmonella_20230928_nie/|while read i; do faprefix=${i//.fa/}; echo -e "sh $PWD/1_mlst_seqsero2.sh ~/data/Salmonella_20230928_nie/${i} ${faprefix} ~/ST_Serotype/ST_Serotype_result">>1_mlst_seqsero2_run.list; done
split -l 100 1_mlst_seqsero2_run.list rundata/1_mlst_seqsero2_run.list.split
ls rundata/1_mlst_seqsero2_run.list.split*|while read i; do echo -e "sh $PWD/${i}">>1_mlst_seqsero2_run.list.split.sh; done
cd slurm_out
sbatch -J st -A p_phage -p batch4 -n 4 -N 1 -a 1-189 ~/script/wrapper.sh ~/ST_Serotype/1_mlst_seqsero2_run.list.split.sh 1
sbatch -J cat -A p_phage -p batch4 -n 20 -N 1 2_catresult_SeqSero2.sh
sbatch -J cat -A p_phage -p batch4 -n 20 -N 1 2_catresult_MLST.sh


############# own ###################################
cd ~/ARG_combine
sh makediamonddb.sh
ls ~/MGE/MGE_result/Salmonella/|while read i; do echo "sh ~/ARG_analyse5_add_mutation_models.sh ~/MGE/MGE_result/Salmonella/${i}/prokka/${i}.faa ~/MGE/MGE_result/Salmonella/${i}/ ${i}">>ARG_analyse5_add_mutation_models_run.list; done
sh ARG_analyse5_add_mutation_models_run.list


############# RefSeq ###############################
cd ~/RefSeq_Anno
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
mv assembly_summary.txt refseq_bacteria_assembly_summary_20240702.txt
awk -F'\t' '{if($8~/Salmonella/){print}}' refseq_bacteria_assembly_summary_20240702.txt>refseq_bacteria_assembly_summary_20240702_Salmonella.txt
awk -F'\t' '{arr=split($20,a,"/"); print "cp",$20"/"a[arr]"_cds_from_genomic.fna.gz ~/RefSeq_Anno/cds_data/"}' refseq_bacteria_assembly_summary_20240702_Salmonella.txt |sed 's#https://ftp.ncbi.nlm.nih.gov/genomes/all#/lustre/public/commondatabase/assembly#g'>1_cp_file.sh
sed -i 's#_protein.faa.gz#_translated_cds.faa.gz#g' 1_cp_file_faa.sh
split -l 50 1_cp_file_faa.sh rundata/1_cp_file_faa.sh.split
ls rundata/1_cp_file_faa.sh.split*|while read i; do echo "sh ${i}">>1_cp_file_faa.sh.split.sh; done
sbatch -J cp -A p_phage -p batch2 -n 1 -N 1 -a 1-285 ~/script/wrapper.sh 1_cp_file_faa.sh.split.sh 1
cat slurm-*|awk -F"’" '{print $1}'|sed 's#cp: cannot stat ‘##g'|sed 's#/lustre/public/commondatabase/assembly#https://ftp.ncbi.nlm.nih.gov/genomes/all#g'|while read i; do echo $i|awk -v i=$i '{arr=split(i,a,"/"); print "wget -O ~/RefSeq_Anno/faa_data/"a[arr]" "i}'; done>>2_wget_file_faa_bc.sh
sh 2_wget_file_faa_bc.sh
perl 3_cat_filt_faa.pl
perl 4_get_changeidinfo_faa.pl
perl 4_2_get_changeidcds_faa.pl
sh 5_0_db_rmdup.sh
sh 5_diamond_makedb.sh

awk '{print "sh 6_diamond_eachgenome.sh ~/MGE/MGE_result/Salmonella/"$1"/prokka/"$1".faa ~/RefSeq_Anno/Anno_Result/"$1".diamond.out"}' ../Salmonella_MIC_usegenome.list>>6_diamond_eachgenome_run.list
split -l 20 6_diamond_eachgenome_run.list rundata/6_diamond_eachgenome_run.list.split
ls rundata/6_diamond_eachgenome_run.list.split*|while read i; do echo "sh ${i}">>6_diamond_eachgenome_run.list.split.sh; done
sbatch -J dia -A p_phage -p batch4 -n 10 -N 1 -a 1-745 ~/script/wrapper.sh 6_diamond_eachgenome_run.list.split.sh 1
sbatch -J combine -A p_phage -p gpu -n 40 -N 1 7_combine_refseq_diamond_result_run.sh
perl 8_add_info_RefSeq_Anno.pl 
perl 9_combine_refseq_prokka_arg_diamond_allfile.pl ~/RefSeq_Anno/RefSeq_Anno_Result_14888file_addinfo_add_prokka_arg_fromdiamond.txt
perl 9_combine_refseq_prokka_arg_diamond_allfile.pl ~/RefSeq_Anno/RefSeq_Anno_Result_14888file_addinfo_add_prokka_arg_fromdiamond_add_mutation_models.txt


##### Generate drug resistance annotation results
#temp0=""
temp0="_add_mutation_models"
awk -F'\t' '{if($25!="NA"){print}}' RefSeq_Anno_Result_14888file_addinfo_add_prokka_arg_fromdiamond${temp0}.txt>RefSeq_Anno_Result_14888file_addinfo_add_prokka_arg_fromdiamond_havearg${temp0}.txt
#temp="RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg"
temp="RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_add_mutation_models"
awk -F'\t' '{if($28>=90 && $29>=80 && $30>=80){print}}' ${temp}.txt>${temp}_PASSarg908080.txt
awk -F'\t' '{if($28<90 || $29<80 || $30<80){print}}' ${temp}.txt>${temp}_notPASSarg908080.txt
awk -F'\t' '{if($12>=90 && $13>=80 && $14>=80){print}}' ${temp}_notPASSarg908080.txt>${temp}Sarg908080_PASSrefseq908080.txt
awk -F'\t' '{if($12<90 || $13<80 || $14<80){print}}' ${temp}_notPASSarg908080.txt>${temp}_notPASSarg908080_notPASSrefseq908080.txt
awk -F'\t' '{print $18"\t"$20}' ${temp}_notPASSarg908080_PASSrefseq908080.txt|sort|uniq>${temp}_notPASSarg908080_PASSrefseq908080_GeneProteinList.txt
awk -F'\t' '{print $20}' ${temp}_notPASSarg908080_PASSrefseq908080.txt|sort|uniq>${temp}_notPASSarg908080_PASSrefseq908080_ProteinList.txt
perl 10_get_final_havearg_notPASSarg908080_PASSrefseq908080.pl ${temp}
awk -F'\t' '{if($12>=40 && $13>=40 && $14>=40){print}}' ${temp}_notPASSarg908080_notPASSrefseq908080.txt>${temp}_notPASSarg908080_notPASSrefseq908080_PASSrefseq404040.txt
awk -F'\t' '{if(NR==1){print $0"\tAnnotation_Source\tUse_Gene_Name\tUse_Protein_Name"}}' ${temp}_PASSarg908080.txt>${temp}_last_use_result.txt
awk -F'\t' 'NR==FNR{str=$4"\t"$6; a[$1]=str}NR>FNR{if($33 in a){print $0"\tARGcombine\t"a[$33]}}' ~/ARG_combine/ARG_combine_info_aa_use2${temp0}.txt ${temp}_PASSarg908080.txt >>${temp}_last_use_result.txt
sed -n '2,$p' ${temp}_notPASSarg908080_PASSrefseq908080_use.txt >>${temp}_last_use_result.txt
awk -F'\t' 'NR==FNR{str=$4"\t"$6; a[$1]=str}NR>FNR{if($33 in a){print $0"\tARGcombine\t"a[$33]}}' ~/ARG_combine/ARG_combine_info_aa_use2${temp0}.txt ${temp}_notPASSarg908080_notPASSrefseq908080_PASSrefseq404040.txt >>${temp}_last_use_result.txt
awk -F'\t' '{print $(NF-1)"\t"$NF"\t"$(NF-2)}' ${temp}_last_use_result.txt|sort -k 1,2|uniq>${temp}_last_use_result_name.txt
perl 13_AllUsedResult_change_GeneProteinName.pl ${temp}
#perl 14_get_ARG_finalResult_14880.pl 1
perl 14_get_ARG_finalResult_14880.pl 2
perl 15_get_MIC_ARG_result.pl ../MIC_20230906.txt 14880file_ARGanno_Protein_result_copynumber${temp0}.txt 14880file_ARGanno_Protein_result_copynumber${temp0}_addPhenotype.txt
perl 15_get_MIC_ARG_result.pl ../MIC_20230906.txt 14880file_ARGanno_Protein_result_absence01${temp0}.txt 14880file_ARGanno_Protein_result_absence01${temp0}_addPhenotype.txt
Rscript 16_get_ARGreference.R ${temp}_last_use_result_last.txt ${temp}_last_use_result_last_selectSNPreference.txt
perl 17_get_ARGSNPreference_fna.pl ${temp}_last_use_result_last_selectSNPreference.txt ${temp}_last_use_result_last_selectSNPreference.ffn


### callSNP
cd ~/ARGSNP
cat ~/Salmonella_MIC_usegenome_20240829.list|while read i
do
        echo "sh ~/ARGSNP/1_mummer_callsnp.sh ~/ARGSNP/mummer${temp0} ~/ARGSNP/${temp}_last_use_result_last_selectSNPreference.ffn ~/MGE/MGE_result/Salmonella/${i}/prokka/${i}.ffn ${i} 10">>1_mummer_callsnp_run.list
done
split -l 200 1_mummer_callsnp_run.list rundata/1_mummer_callsnp_run.list.split
ls rundata/1_mummer_callsnp_run.list.split*|while read i; do echo "sh $PWD/${i}">>1_mummer_callsnp_run.list.split.sh; done
sbatch -J callsnp -A p_phage -p batch4 -n 10 -N 1 -a 1-75 ~/script/wrapper.sh 1_mummer_callsnp_run.list.split.sh 1
cat ../Salmonella_MIC_usegenome_20240829.list |while read i
do
	echo -e "perl $PWD/2_select_realARGgene_SNPresult.pl ~/ARGSNP/mummer${temp0}/${i}.snps.snp.indel959090 ${i} ~/ARGSNP/mummer${temp0}/${i}.snps.snp.indel959090.filtedARG">>2_select_realARGgene_SNPresult_run.list
done
split -l 500 2_select_realARGgene_SNPresult_run.list rundata/2_select_realARGgene_SNPresult_run.list.split
ls rundata/2_select_realARGgene_SNPresult_run.list.split*|while read i; do echo "sh $PWD/${i}">>2_select_realARGgene_SNPresult_run.list.split.sh; done
sbatch -J callsnp -A p_phage -p batch3 -x node26 -n 1 -N 1 -a 1-90 ~/script/wrapper.sh 2_select_realARGgene_SNPresult_run.list.split.sh 1
##id95cov90---01
perl 3_get_14880_SNPmatrix.pl ~/RefSeq_Anno/${temp}_last_use_result_last_selectSNPreference.txt ~/ARGSNP/mummer${temp0} ${temp}_last_use_result_last_SNPmatrix959090.txt
perl 4_get_MIC_SNP_result.pl ../MIC_20230906.txt ${temp}_last_use_result_last_SNPmatrix959090.txt ${temp}_last_use_result_last_SNPmatrix959090_addPhenotype.txt
##id95cov90---YesNoNA
perl 3_2_get_14880_SNPmatrix_YesNoNA.pl ~/RefSeq_Anno/${temp}_last_use_result_last_selectSNPreference.txt ~/ARGSNP/mummer${temp0} ${temp}_last_use_result_last_SNPmatrix959090_YesNoNA.txt
perl 4_get_MIC_SNP_result.pl ../MIC_20230906.txt ${temp}_last_use_result_last_SNPmatrix959090_YesNoNA.txt ${temp}_last_use_result_last_SNPmatrix959090_YesNoNA_addPhenotype.txt

######## ARG combine SNP ##########
cd ~/ARG_with_SNP
Rscript 1_combine_arg_snpindel.R 

######## MGE anno ######
ls ~/data/Salmonella_20230928_nie/|while read i; do faprefix=${i//.fa/}; echo -e "sh ~/Feature_analysis/MGE_gcPathogen_bacteria_pipeline.sh -O ~/MGE_result/ -t 20 -I ~/data/Salmonella_20230928_nie/${i} -i 90 -s 80 -q 80 -e 0.00001">>MGE_gcPathogen_bacteria_pipeline_run.list; done
split -l 100 MGE_gcPathogen_bacteria_pipeline_run.list rundata/MGE_gcPathogen_bacteria_pipeline_run.list.split
ls rundata/MGE_gcPathogen_bacteria_pipeline_run.list.split*|while read i;do echo -e "sh $PWD/$i">>MGE_gcPathogen_bacteria_pipeline_run.list.split.sh;done
sbatch -J mge -A p_phage -p batch3 -n 20 -N 1 -a 1-149 ~/script/wrapper.sh MGE_gcPathogen_bacteria_pipeline_run.list.split.sh 1


############## plublic genomes ###############
awk -F'\t' '{print "sh ~/Feature_analysis/1_all_ARGcombine_RefSeq_anno.sh "$3"/prokka/"$2".faa ~/Feature_analysis/ARGcombine_RefSeq_result/ "$1}' Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort.txt>Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list
split -l 50 Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list rundata/Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split
ls rundata/Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split*|while read i
do
        echo "sh $PWD/${i}">>Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh
done
sed -n '1,1000p' Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh >Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_1_1000
sed -n '1001,2000p' Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh >Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_1001_2000
sed -n '2001,3000p' Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh >Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_2001_3000
sed -n '3001,$p' Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh >Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_3001_3686
sbatch -J abricat -A p_phage -p batch4 -n 10 -N 1 -a 1-1000 ~/script/wrapper.sh ../Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_1_1000 1
sbatch -J abricat -A p_phage -p batch2 -n 10 -N 1 -a 1-1000 ~/script/wrapper.sh ../Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_1001_2000 1
sbatch -J abricat -A p_phage -p batch2 -n 10 -N 1 -a 1-1000 ~/script/wrapper.sh ../Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_2001_3000 1
sbatch -J abricat -A p_phage -p batch2 -n 10 -N 1 -a 1-686 ~/script/wrapper.sh ../Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_runARGanno.list.split.sh_3001_3686 1
##same as aforementioned.

####### plasmid mobility predict use MOB_suite.
cd ~/MGE
ls ~/data/Salmonella_20230928_nie|sed 's#.fa$##g'|while read i;do echo -e "sh $PWD/1_get_plasmid.sh $i">>1_get_plasmid_run.list;done
split -l 100 1_get_plasmid_run.list rundata/1_get_plasmid_run.list.split
ls rundata/1_get_plasmid_run.list.split*|while read i;do echo -e "sh $PWD/$i">>1_get_plasmid_run.list.split.sh;done
sbatch -J a -A p_phage -p batch4 -n 4 -N 1 -a 1-149 ~/script/wrapper.sh 1_get_plasmid_run.list.split.sh 1
ls ~/data/Salmonella_20230928_nie|sed 's#.fa$##g'|while read i;do echo -e "sh $PWD/2_MOB_suite.sh $i 8">>2_MOB_suite_run.list;done
split -l 100 2_MOB_suite_run.list rundata/2_MOB_suite_run.list.split
ls rundata/2_MOB_suite_run.list.split*|while read i;do echo -e "sh $PWD/$i">>2_MOB_suite_run.list.split.sh;done
sbatch -J a -A p_phage -p batch4 -n 4 -N 1 -a 1-149 ~/script/wrapper.sh 2_MOB_suite_run.list.split.sh 1
#combine 14880 mobsuite result in ~/MGE/MGE_result/Salmonella_MOBsuite/14880_MOB_suite_result.txt
perl 3_MOB_suite_result_merge.pl


###################################### Model ############################################
cd ~/ARG_MIC
sbatch -J combine -A p_phage -p batch4 -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 0_split_RSgroup_82_run.list 1
sbatch -J reauc -A p_phage -p batch2 -n 20 -N 1 -a 1-13 ~/script/wrapper.sh 1_sklearn_8cv10fold_default_allfeature_absence01_verifytest_combine_run.list 1
sbatch -J reauc -A p_phage -p batch2 -n 20 -N 1 -a 1-13 ~/script/wrapper.sh 1_sklearn_8cv10fold_default_allfeature_copynumber_verifytest_combine_run.list 1
sbatch -J all10absence -A p_phage -p batch2 -n 20 -N 1 -a 1-13 ~/script/wrapper.sh 2_sklearn_8cv10fold_default_allfeature_absence01_topfeature_verifytest_combine_run.list 1
sbatch -J all10absence -A p_phage -p batch2 -n 20 -N 1 -a 1-14 ~/script/wrapper.sh 2_sklearn_8cv10fold_default_allfeature_copynumber_topfeature_verifytest_combine_run.list 1

cd ~/ARGSNP
sbatch -J combine -A p_phage -p batch3 -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 5_split_RSgroup_82_run.list 1
sbatch -J combine -A p_phage -p batch3 -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 5_split_RSgroup_82_YesNoNA_run.list 1
sbatch -J combine -A p_phage -p batch2 -n 20 -N 1 -a 1-13 ~/script/wrapper.sh 6_sklearn_8cv10fold_default_allfeature_combine_verifytest_run.list 1
sbatch -J retrain -A p_phage -p batch2 -n 20 -N 1 -a 1-13 ~/script/wrapper.sh 6_sklearn_8cv10fold_default_allfeature_YesNoNA_categoryYesNoNA_and_combine_verifytest_run.list 1

cd ~/ARG_with_SNP
Rscript 1_combine_arg_snpindel.R
sbatch -J combine -A p_phage -p batch2 -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 2_split_RSgroup_811_ARGanno_absence01_SNPmatrix959090_YesNoNA_run.list 1
sbatch -J combine -A p_phage -p batch3 -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 2_split_RSgroup_811_ARGanno_absence01_SNPmatrix959090_absence01_run.list 1
sbatch -J combine -A p_phage -p batch3 -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 2_split_RSgroup_811_ARGanno_copynumber_SNPmatrix959090_YesNoNA_run.list 1
sbatch -J combine -A p_phage -p batch3 -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 2_split_RSgroup_811_ARGanno_copynumber_SNPmatrix959090_absence01_run.list 1
sbatch -J retrain -A p_phage -p batch4 -n 40 -N 1 -a 1-13 ~/script/wrapper.sh 3_sklearn_8cv10fold_default_allfeature_ARGanno_absence01_SNPmatrix959090_absence01_category_and_combine_verifytest_run.list 1
sbatch -J retrain -A p_phage -p batch4 -n 40 -N 1 -a 1-13 ~/script/wrapper.sh 3_sklearn_8cv10fold_default_allfeature_ARGanno_copynumber_SNPmatrix959090_absence01_category_and_combine_verifytest_run.list 1
sbatch -J retrain -A p_phage -p batch4 -n 40 -N 1 -a 1-13 ~/script/wrapper.sh 3_sklearn_8cv10fold_default_allfeature_ARGanno_absence01_SNPmatrix959090_YesNoNA_category_and_combine_verifytest_run.list 1
sbatch -J retrain -A p_phage -p batch4 -n 40 -N 1 -a 1-13 ~/script/wrapper.sh 3_sklearn_8cv10fold_default_allfeature_ARGanno_copynumber_SNPmatrix959090_YesNoNA_category_and_combine_verifytest_run.list 1
sbatch -J topfeature -A p_phage -p batch2 -n 40 -N 1 -a 1-7 ~/script/wrapper.sh 4_sklearn_8cv10fold_default_allfeature_ARGanno_absence01_SNPmatrix959090_absence01_topfeature_category_and_combine_verifytest_run.list 1
sbatch -J topfeature -A p_phage -p batch2 -n 40 -N 1 -a 1-6 ~/script/wrapper.sh 4_sklearn_8cv10fold_default_allfeature_ARGanno_absence01_SNPmatrix959090_YesNoNA_topfeature_category_and_combine_verifytest_run.list 1
sbatch -J topfeature -A p_phage -p batch2 -n 40 -N 1 -a 1 ~/script/wrapper.sh 4_sklearn_8cv10fold_default_allfeature_ARGanno_absence01_SNPmatrix959090_absence01_topfeature_category_and_combine_verifytest_nooptimal_run.list 1
sbatch -J topfeature -A p_phage -p batch2 -n 40 -N 1 -a 1-2 ~/script/wrapper.sh 4_sklearn_8cv10fold_default_allfeature_ARGanno_absence01_SNPmatrix959090_YesNoNA_topfeature_category_and_combine_verifytest_nooptimal_run.list 1
sbatch -J topfeature -A p_phage -p batch2 -n 20 -N 1 -a 1-4 ~/script/wrapper.sh 4_sklearn_8cv10fold_default_allfeature_ARGanno_copynumber_SNPmatrix959090_absence01_topfeature_category_and_combine_verifytest_run.list 1
sbatch -J topfeature -A p_phage -p batch2 -n 20 -N 1 -a 1 ~/script/wrapper.sh 4_sklearn_8cv10fold_default_allfeature_ARGanno_copynumber_SNPmatrix959090_YesNoNA_topfeature_category_and_combine_verifytest_run.list 1




################### Use optimalized hyperparameter, to train topfeature models in all training set, and get final models of antibiotic, as well as feature importance and evaluation results.
cd ~/Feature_analysis
python 12_XGBoost_topfeature_final_module_AMP_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_AMS_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_AZM_ARGMIC_absence01_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_CAZ_ARGMIC_absence01_SNP_absence01.py
python 12_XGBoost_topfeature_final_module_CFZ_ARGMIC_copynumber_SNP_absence01.py
python 12_XGBoost_topfeature_final_module_CHL_ARGMIC_copynumber.py
python 12_XGBoost_topfeature_final_module_CIP_ARGMIC_absence01_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_CT_ARGMIC_absence01_SNP_absence01.py
python 12_XGBoost_topfeature_final_module_CTX_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_GEN_ARGMIC_absence01_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_NAL_ARGMIC_absence01_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_SXT_ARGMIC_absence01_SNP_absence01.py
python 12_XGBoost_topfeature_final_module_TET_ARGMIC_absence01.py
#### final models is ~/Feature_analysis/use_topfeature/XGBoost*_topfeature_finalmodel_traindata.txt, however, when using it specifically, it is necessary to refer to the processing and training methods for this training set matrix in the 13_XGBoost_topfeature_final_madule_XX_predict.py file to train the model and predict.


####### auc plot
python ~/Plot/3_plot_roc_add_mutaion_models.py


###################### Global risk analysis ###################################################################################################
cd ~/Feature_analysis
for i in AMP AMS AZM CIP CT CTX CAZ SXT NAL GEN CIP SXT TET
do
	python 13_XGBoost_topfeature_final_module_${i}_predict.py
done

## Add mobility
awk -F'\t' 'NR==FNR{arrnum=split($2,arr," ");str=$1"\t"arr[1];a[str]=$15}NR>FNR{if($4~/plasmid/){str2=$2"\t"$5;if((str2 in a) && (a[str2]!~/non-mobilizable/)){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"a[str2]}}else{print}}' ~/MGE/MOB_suite_result.txt ~/Feature_analysis/updownstreamInfo_202605/all_species_MGEresult_filtANI96.txt>~/Feature_analysis/updownstreamInfo_202605/all_species_MGEresult_filtANI96_use0.txt
awk -F'\t' 'NR==FNR{arrnum=split($2,arr," ");str=$1"\t"arr[1];a[str]=$7}NR>FNR{if($4~/plasmid/){str2=$2"\t"$5;if((str2 in a) && $11=="-"){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"a[str2]"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18}else{print }}else{print}}' ~/MGE/MOB_suite_result.txt ~/Feature_analysis/updownstreamInfo_202605/all_species_MGEresult_filtANI96_use0.txt>~/Feature_analysis/updownstreamInfo_202605/all_species_MGEresult_filtANI96_use.txt


######### re statistic MGEs's fareast lineage --add plasmid mobility
cd ~/Feature_analysis
awk -F'\t' -v OFS="\t" 'NR==FNR{print $0;a[$1]=1}NR>FNR{if($1!="genbank-accession"){if(!($1 in a)){print $0"\tSalmonella_enterica"}}}' ~/MGE/all_species_runned_genome_ANI96need_genomelist.txt Salmonella_public_184251_ANI96need_genomelist.txt>All_MGE_use_public_ANI96need_genomelist.txt
perl 15_get_eachMGE_species_prop.pl All_MGE_use_public_ANI96need_genomelist.txt updownstreamInfo_202605/all_species_MGEresult_filtANI96_use.txt updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_crossSpecies_prop.txt
perl 16_add_feature_mgeallspecies_lineage.pl updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_crossSpecies_prop.txt ~/MGE/bact.lineage2.havemge.csv updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_crossSpecies_prop_addlineage.txt
perl 17_count_eachMGE_crossedLineage_genome2.pl ~/MGE/all_species_runned_lineagenumber.txt updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_crossSpecies_prop_addlineage.txt updownstreamInfo_202605/all_species_MGEresult_filtANI96_use_allMGE_crossedLineage_genome1.txt
perl 17_2_count_eachMGE_crossedLineage_selectuse_genome1use_202605.pl
#removed plasmid -       Phylum  Proteobacteria;Firmicutes;Bacteroidetes;Actinobacteria;Actinomycetota;Fusobacteria, which will will mislead the calculation resultswill mislead the calculation results

#####Calculate the correlation between features and MGE in the genome
## own genomes mge asso
cd ~/Feature_analysis
cat ../Salmonella_MIC_usegenome_20240829.list |while read i; do echo -e "perl $PWD/23_get_own14880_argAssomge_eachgenome.pl ~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_add_mutation_models_last_use_result_last.txt ~/MGE/MGE_result/Salmonella/${i}/${i}_MGE_merged_site.txt ~/MGE/MGE_result/Salmonella/${i}/${i}_MGE_merged_ARGcombine_result_add_mutation_models.txt";echo -e "perl $PWD/23_get_own14880_argAssomge_eachgenome.pl ~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last.txt ~/MGE/MGE_result/Salmonella/${i}/${i}_MGE_merged_site.txt ~/MGE/MGE_result/Salmonella/${i}/${i}_MGE_merged_ARGcombine_result.txt">>23_get_own14880_argAssomge_eachgenome_run.list;done
split -l 100 23_get_own14880_argAssomge_eachgenome_run.list rundata/23_get_own14880_argAssomge_eachgenome_run.list.split
ls rundata/23_get_own14880_argAssomge_eachgenome_run.list.split*|while read i
do
        echo "sh $PWD/${i}">>23_get_own14880_argAssomge_eachgenome_run.list.split.sh
done
sbatch -J mgemerge -A p_phage -p batch24 -n 4 -N 1 -a 1-298 ~/script/wrapper.sh 23_get_own14880_argAssomge_eachgenome_run.list.split.sh 1

## public genomes mge asso
cat MGE_Salmonella_enterica_runedgenome_path_2014use.txt|while read i; do arr=($i); echo -e "perl $PWD/18_get_argAssomge_eachgenome.pl /lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/ARGcombine_RefSeq_result/${arr[0]}/${arr[0]}_Salmonella_RefSeq_add_prokka_arg_fromdiamond_last_use_result.txt ${arr[1]}/${arr[0]}_MGE_merged_site.txt $PWD/MGE_prerunned/${arr[0]}_MGE_merged_ARGcombine_result.txt";echo -e "perl $PWD/18_get_argAssomge_eachgenome.pl /lustre/home/niejingyi2023/project/bjCDC/Feature_analysis/ARGcombine_RefSeq_result/${arr[0]}/${arr[0]}_Salmonella_RefSeq_add_prokka_arg_fromdiamond_add_mutation_models_last_use_result.txt ${arr[1]}/${arr[0]}_MGE_merged_site.txt $PWD/MGE_prerunned/${arr[0]}_MGE_merged_ARGcombine_result_add_mutation_models.txt">>18_get_public184251_argAssomge_eachgenome_run.list;done
split -l 400 18_get_public184251_argAssomge_eachgenome_run.list rundata/18_get_public184251_argAssomge_eachgenome_run.list.split
ls rundata/18_get_public184251_argAssomge_eachgenome_run.list.split*|while read i
do
        echo "sh $PWD/${i}">>18_get_public184251_argAssomge_eachgenome_run.list.split.sh
done
sbatch -J mgemerge -A p_phage -p batch4 -x node44 -n 2 -N 1 -a 1-922 ~/script/wrapper.sh 18_get_public184251_argAssomge_eachgenome_run.list.split.sh 1

## merge genome feature and mge position
sbatch -J mgemerge -A p_phage -p batch4 -n 40 -N 1 -a 1-8 ~/script/wrapper.sh 19_merge_mgeARGinfo_addpos_202605_run.sh 1

###### genome associated with MGE's information
sbatch -J a -A p_phage -p batch -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 29_2use_get_genome_usemge_info_202605_run.list 1
sbatch -J a -A p_phage -p batch -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 32_get_Rgenome_202605_run.list 1

########## risk score
sbatch -J a -A p_phage -p batch -n 10 -N 1 -a 1-13 ~/script/wrapper.sh 33_risk_score_202605_run.list 1
#### combine antibiotics score
perl 33_2_combine_risk_score_202605.pl

### risk level
awk -F'\t' '{if($1=="Genome"){printf $0}else{printf "\n"$1;for(i=2;i<=NF;i++){if($i=="-"){printf "\t-"}else if($i>=12){printf "\tI"}else if($i>=4){printf "\tII"}else{printf "\tIII"}}}}' predict_use_data_add_mutation_models/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_202605.txt>predict_use_data_add_mutation_models/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel_202605.txt
awk -F'\t' 'NR==FNR{a[$1]=$(NF-2)"\t"$NF}NR>FNR{if($1 in a){print a[$1]"\t"$0}}' predict_use_data_add_mutation_models/Salmonella_enterica_ownuse14846_publicANI96_173348_GEN_topfeature_matrix_add_year_host_country.txt predict_use_data_add_mutation_models/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel_202605.txt>predict_use_data_add_mutation_models/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel_add_year_country_202605.txt




