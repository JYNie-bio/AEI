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
python ~/Plot/26_plot_roc_add_mutaion_models.py



###################### Global risk analysis ###################################################################################################
cd ~/Feature_analysis
for i in AMP AMS AZM CIP CT CTX CAZ SXT NAL GEN CIP SXT TET
do
	python 13_XGBoost_topfeature_final_module_${i}_predict.py
done





##### 整合耐药基因结果生成拷贝数和有无矩阵 ###
sbatch -J matrix -A p_phage -p gpu -n 40 -N 1 10_get_ARG_protein_copynumber_absence01_matrix_run.sh 
#4902325
##### 整合snp结果生成有无矩阵 ###
因为topfeature结果中没有用SNP的YesNoNA的，都是01，因此本处只处理了01的snp结果
sbatch -J matrix -A p_phage -p batch3 -n 40 -N 1 10_get_ARG_SNP_absence01_matrix_run.sh
#4902326
##### 整合snp结果生成YesNoNA矩阵 ###
perl 10_get_ARG_SNP_YesNoNA_matrix.pl 



##### 生成每种抗生素184251个基因组top特征的模型输入矩阵
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R TET absence01 nosnp
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R AZM absence01 snpyesnona
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CAZ absence01 snpabsence
#~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CFX absence01 snpabsence
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CFX absence01 nosnp
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CT copynumber snpyesnona
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CTX absence01 nosnp
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R NAL copynumber snpyesnona
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R AMP absence01 nosnp
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R AMS absence01 nosnp
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CFZ copynumber snpabsence
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CHL copynumber nosnp
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R CIP absence01 nosnp
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R GEN absence01 snpyesnona
~/miniconda3/envs/r/bin/Rscript 11_select_topfeature_184251_matrix.R SXT absence01 snpabsence

注意原校正后feature名字中cryptic_aminoglycoside_nucleotidyltransferase_ANT(3'')/ANT(9)更名为ANT(3'')-Ia_family_aminoglycoside_nucleotidyltransferase_AadA，不影响结果
#注意原校正后feature名字中aminoglycoside_6'-N-acetyltransferase_AAC(6')-Iy更名为aminoglycoside_N-acetyltransferase_AAC(6')-Iy，不影响结果
#注意feature名字中tetracyline_resistance-associated_transcriptional_repressor_TetC和tetracycline_efflux_MFS_transporter_Tet(C)当作不是同一个基因对待的，一个是转录抑制因子，一个是转运蛋白，如果更改的话，会影响涉及相关SNP的结果（如NAL抗生素）。因为原来14880基因组中注释到tetracycline_efflux_MFS_transporter_Tet(C)的只有一条。
### 有一些名字不同的没有更改，因为不在topfeature中
### 原校正后feature名字中aminoglycoside_6'-N-acetyltransferase_AAC(6')-Iy更名为aminoglycoside_N-acetyltransferase_AAC(6')-Iy
### awk -F'\t' 'NR==FNR{str=$1"||"$2"||"$3;a[str]=$1"||"$2"||"$3"||"$4"||"$5}NR>FNR{str2=$1"||"$2"||"$3;if(str2 in a){if(a[str2]!=$1"||"$2"||"$3"||"$4"||"$5){print $0}}}' RefSeq_Anno/RefSeq_Anno_Result_14888file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_needChangeInfo.txt Feature_analysis/Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_arg_last_use_result_needChangeInfo.txt

检查11_select_topfeature_184251_matrix.R结果
###检查后没有问题，有两种topfeature中SNP没有在公共数据中注释到，都用0补齐了，见CAZ和NAL抗生素。
#head -1 predict_use_data/Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_arg_last_use_result_last_NAL_topfeature_matrix.txt |awk -F'\t' '{for(i=2;i<=NF;i++){print $i}}'>tempa.txt



##### 生成每种抗生素最终的模型
#Rscript 12_提取top特征的训练集测试集.R 
python 12_XGBoost_topfeature_final_module_AMP_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_AMS_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_AZM_ARGMIC_absence01_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_CAZ_ARGMIC_absence01_SNP_absence01.py
#python 12_XGBoost_topfeature_final_module_CFX_ARGMIC_absence01_SNP_absence01.py
python 12_XGBoost_topfeature_final_module_CFX_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_CFZ_ARGMIC_copynumber_SNP_absence01.py
python 12_XGBoost_topfeature_final_module_CHL_ARGMIC_copynumber.py
python 12_XGBoost_topfeature_final_module_CIP_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_CT_ARGMIC_copynumber_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_CTX_ARGMIC_absence01.py
python 12_XGBoost_topfeature_final_module_GEN_ARGMIC_absence01_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_NAL_ARGMIC_copynumber_SNP_YesNoNA.py
python 12_XGBoost_topfeature_final_module_SXT_ARGMIC_absence01_SNP_absence01.py
python 12_XGBoost_topfeature_final_module_TET_ARGMIC_absence01.py

##### 整合每种抗生素最终模型的效果
head -1 use_topfeature/XGBoost_GEN_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'\t' '{print "Antibiotic\tML_data\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'>use_topfeature/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt
ls use_topfeature/XGBoost_*_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'_' '{print $3}'|sort|uniq|while read i; do sed -n '2,$p' use_topfeature/XGBoost_${i}_topfeature_finalmodel_evaluation_indicator_results.txt|awk -v i=$i '{print i"\t"$0}'>>use_topfeature/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt; done

###发现抗生素CT用最优超参数计算后，测试集和独立测试集效果不好，怀疑最优超参数只在训练集上训练寻找可能找到的不是最优的模型，因此比较默认超参数的模型和最优超参数模型，差不多的情况下选择最优超参数模型。
head -1 use_topfeature2/XGBoost_GEN_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'\t' '{print "Antibiotic\tML_data\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'>use_topfeature2/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt
ls use_topfeature2/XGBoost_*_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'_' '{print $3}'|sort|uniq|while read i; do sed -n '2,$p' use_topfeature2/XGBoost_${i}_topfeature_finalmodel_evaluation_indicator_results.txt|awk -v i=$i '{print i"\t"$0}'>>use_topfeature2/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt; done


sbatch -J test -A p_phage -p batch4 -n 10 -N 1 -a 1-14 ~/script/wrapper.sh 13_test_run.list 1
#5104859
head -1 use_topfeature_test/XGBoost_GEN_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'\t' '{print "Antibiotic\tML_data\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'>use_topfeature_test/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt
ls use_topfeature_test/XGBoost_*_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'_' '{print $4}'|sort|uniq|while read i; do sed -n '2,$p' use_topfeature_test/XGBoost_${i}_topfeature_finalmodel_evaluation_indicator_results.txt|awk -v i=$i '{print i"\t"$0}'>>use_topfeature_test/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt; done
head -1 use_topfeature_test2/XGBoost_GEN_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'\t' '{print "Antibiotic\tML_data\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'>use_topfeature_test2/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt
ls use_topfeature_test2/XGBoost_*_topfeature_finalmodel_evaluation_indicator_results.txt|awk -F'_' '{print $4}'|sort|uniq|while read i; do sed -n '2,$p' use_topfeature_test2/XGBoost_${i}_topfeature_finalmodel_evaluation_indicator_results.txt|awk -v i=$i '{print i"\t"$0}'>>use_topfeature_test2/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt; done
#diff use_topfeature/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt use_topfeature_test/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt|less -S
#diff use_topfeature2/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt use_topfeature_test2/XGBoost_topfeature_finalmodel_evaluation_indicator_results_Allantibiotic.txt|less -S
#confirmed the model traindata and testdata can be used directly.

sbatch -J roc -A p_phage -p gpu -n 20 -N 1 26_plot_roc_run.sh 
#5105296

##### 预测184251个基因组各抗生素的耐药表型
python 13_XGBoost_topfeature_final_module_TET_predict.py
python 13_XGBoost_topfeature_final_module_AMP_predict.py
python 13_XGBoost_topfeature_final_module_AMS_predict.py
python 13_XGBoost_topfeature_final_module_AZM_predict.py
python 13_XGBoost_topfeature_final_module_CAZ_predict.py
python 13_XGBoost_topfeature_final_module_CFX_predict.py
python 13_XGBoost_topfeature_final_module_CFZ_predict.py
python 13_XGBoost_topfeature_final_module_CHL_predict.py
python 13_XGBoost_topfeature_final_module_CIP_predict.py
python 13_XGBoost_topfeature_final_module_CT_predict.py
python 13_XGBoost_topfeature_final_module_CTX_predict.py
python 13_XGBoost_topfeature_final_module_GEN_predict.py
python 13_XGBoost_topfeature_final_module_NAL_predict.py
python 13_XGBoost_topfeature_final_module_SXT_predict.py
sbatch -J predict -A p_phage -p batch3 -n 10 -N 1 -a 1-14 ~/script/wrapper.sh 13_XGBoost_topfeature_final_module_run.list 1
#5113158


##### 合并各抗生素模型各特征及其重要性
echo -e "Antibiotic\tFeature\tImportance">use_topfeature/XGBoost_topfeature_finalmodel_feature_importances_Allantibiotic.txt
ls use_topfeature/XGBoost_*_topfeature_finalmodel_feature_importances.txt|awk -F'_' '{print $3}'|while read i;do `awk -v i=$i -F'\t' '{if($2!="Feature"){print i"\t"$2"\t"$3}}' use_topfeature/XGBoost_${i}_topfeature_finalmodel_feature_importances.txt>>use_topfeature/XGBoost_topfeature_finalmodel_feature_importances_Allantibiotic.txt`;done



########################################################## 各个模型特征的上下游基因情况--各类MGE跨种属情况 #####################################################
### 各个基因组关联cds区和MGE结果
############################ MGE融合cds #############################
cd ~
ls MGE/MGE_result/Salmonella|while read i; do echo -e "perl ~/script/get_MGEmerged_eachsoftware_prokkacds13_salmonella_use_ISonlyonesite.pl ~/MGE/MGE_result/Salmonella/${i}/prokka/${i}.gff ~/MGE/MGE_result/Salmonella/${i}/ ${i}_MGE_merged_site.txt ${i}_MGE_merged_cds_result_ISonlyonesite.txt">>MGE_mergecds_command_ISonlyonesite.list; done
split -l 100 MGE_mergecds_command_ISonlyonesite.list MGE_mergecds_command_ISonlyonesite.list.split
ls MGE_mergecds_command_ISonlyonesite.list.split*|while read i; do echo "sh $i">>MGE_mergecds_command_ISonlyonesite.list.split.sh; done
sbatch -J mergecds -A p_phage -p gpu -n 2 -N 1 -a 1-169 ~/script/wrapper.sh MGE_mergecds_command_ISonlyonesite.list.split.sh 1
#4853804





##################### 各类MGE跨种属情况 ######################
cd ~/MGE
#合并所有病原菌mge注释结果
perl get_all_species_MGEresult.pl
#生成所有病原菌mge的物种列表
perl get_all_species.pl
#生成所有病原菌跑过mge的物种列表
perl get_all_species_MGErungenome.pl
#生成所有病原菌跑过mge的种所属的各级谱系情况，每级谱系的数目，以及所属此谱系的次级谱系数目，比如跑过mge的属的数目，每种属下种的数目。
perl get_all_species_MGErunned_lineagenumber.pl
#生成目前仍没有mge注释的所用2014年后沙门公共基因组列表
grep "Salmonella_enterica" ~/project/MGE/all_species_runned_genome.list |cut -f2|sort|uniq>~/project/MGE/MGE_Salmonella_enterica_runedgenome.list
awk -F'\t' 'NR==FNR{a[$1]=1}NR>FNR{if($1 in a){print $0}}' ~/project/MGE/MGE_Salmonella_enterica_runedgenome.list Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort.txt >Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_nowhavemge.txt
awk -F'\t' 'NR==FNR{a[$1]=1}NR>FNR{if(!($1 in a)){print $0}}' ~/project/MGE/MGE_Salmonella_enterica_runedgenome.list Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort.txt >Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_needrunmge.txt


cd  ~
#NCBI lineage
awk -F',' '{print $NF}' ~/project/bjCDC/MGE/bact.lineage.csv |sed 's# #_#g'|sed 's#^M##g'>~/project/bjCDC/MGE/species.list
awk 'NR==FNR{a[$1]=$1}NR>FNR{if(!($1 in a)){print $1}}' MGE/species.list ~/project/MGE/all_species_MGEanno_species.list|head
#发现有两个物种命名~/project/bjCDC/MGE/bact.lineage.csv和~/project/MGE/all_species_MGEanno_species.list不一致，Clostridium_symbiosum应命名为[Clostridium]_symbiosum，Haemophilus_ducreyi应命名为[Haemophilus]_ducreyi，但为了统计方便
，将~/project/bjCDC/MGE/bact.lineage.csv中[Clostridium]_symbiosum改为Clostridium_symbiosum，[Haemophilus]_ducreyi改为Haemophilus_ducreyi，写入~/project/bjCDC/MGE/bact.lineage2.csv
cp ~/project/bjCDC/MGE/bact.lineage.csv ~/project/bjCDC/MGE/bact.lineage2.csv
sed -i 's#\[Clostridium\] symbiosum#Clostridium symbiosum#g' ~/project/bjCDC/MGE/bact.lineage2.csv
sed -i 's#\[Haemophilus\] ducreyi#Haemophilus ducreyi#g' ~/project/bjCDC/MGE/bact.lineage2.csv
#提取mge注释的物种谱系
awk -F',' -v OFS="," 'NR==FNR{a[$1]=$1}NR>FNR{gsub(" ","_",$9);if($9 in a){print $0}}' ~/project/MGE/all_species_MGEanno_species.list ~/project/bjCDC/MGE/bact.lineage2.csv>~/project/bjCDC/MGE/bact.lineage2.havemge.csv


#####################################################################################################################################################################################
################################################################################# 后面有加上质粒重新运行，此处结果移动为updownstreamInfo_pre20241204 #################################
############### 统计MGEs的最高可跨种属。##################################
cd ~/Feature_analysis
#统计所有mge的转移情况，不区分时间，mge只统计有具体名字的，IS/integron/transposon/ICE，使用软件MobileElementFinder,DANMEL-blast,BacAnt.
#之后对于每一个沙门基因组，计算ARG关联的MGEs，从上面各个mge的转移情况提取结果即可。

#提取所有病原菌注释MobileElementFinder,DANMEL-blast,BacAnt软件结果
head -1 ~/project/MGE/all_species_MGEresult.txt >updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt.txt
awk -F'\t' -v OFS="\t" '{if($3=="MobileElementFinder"){if($4~/transposon/){str="transposon"}else{str=$4}if($4=="IS"){print $1,$2,$3,str,$5,$6,$7,$8,$9,$10,$12,$11,$13,$14,$15,$16,$17,$18}else{print $1,$2,$3,str,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}}}' ~/project/MGE/all_species_MGEresult.txt >>updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt.txt
awk -F'\t' '{if($3=="DANMEL" || $3=="BacAnt"){print}}' ~/project/MGE/all_species_MGEresult.txt >>updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt.txt
#生成所有mge所在物种。添加某物种某mge的基因组数目，以及此物种所有基因组数目。
perl 15_get_eachMGE_species_prop.pl ~/MGE/all_species_runned_genome.list updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt.txt updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_crossSpecies_prop.txt
#添加七级lineage
perl 16_add_feature_mgeallspecies_lineage.pl updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_crossSpecies_prop.txt ~/MGE/bact.lineage2.havemge.csv updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_crossSpecies_prop_addlineage.txt
#生成所有mge在所有病原菌注释中可跨越谱系，生成七级lineage各级运行mge的数目，比如某科下运行mge的属的种类数目（可跨越的谱系如果只有1个，但不是Salmonella_enterica所属的，也算可跨越谱系）
###只有MGE出现在物种中至少2%（找不到相关文献支持，2%，5%, 10%的都定义看看结果吧）的基因组中才算可移动至此物种。（发现2%时移动mge（可跨谱系传播）占所有mge约35%，5%时占12.5%，因此选用2%可能更合理）
perl 17_count_eachMGE_crossedLineage_num.pl ~/MGE/all_species_runned_lineagenumber.txt updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_crossSpecies_prop_addlineage.txt updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02.txt 0.02
#生成所有mge在所有病原菌注释中最高可跨越谱系等级
#比如定义跨属，只要至少跨两个属就可以，只要属中有一个物种出现即可。
perl 17_2_count_eachMGE_crossedLineage_selectuse.pl 
#生成是否在WHO优先病原体清单里的MGE结果
echo -e "MGE_Type\tElement_Type\tCrossedLineage\tCrossedName\tCrossedLineageNumber\tLineageAllNumber\tCrossedLineageProportion\tWHOpriority">updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority.txt
awk -F'\t' 'NR==FNR{if($1!="Pathogen_name"){a[$1]=$3}}NR>FNR{arr=split($4,arr2,/;/);for(i=1;i<=arr;i++){if(arr2[i] in a){print $0"\t"a[arr2[i]]}}}' WHO_priority_pathogen_list.txt updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02.txt>>updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority0.txt
#合并WHO结果和跨谱系结果
##echo -e "MGE_Type\tElement_Type\tCrossedLineage\tCrossedName\tWHOpriority">updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority1.txt
##awk -F'\t' '{if($1!="MGE_Type"){str=$1"\t"$2"\t"$3"\t"$4;if(!(str in a)){a[str]=$8}else{a[str]=a[str]";"$8}}}END{for(str2 in a){if(a[str2]~/Critical priority/){print str2"\tCritical priority"}else if(a[str2]~/High priority/){print str2"\tHigh priority"}else if(a[str2]~/Medium priority/){print str2"\tMedium priority"}else{print str2"\t-"}}}' updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority0.txt>>updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority.txt
echo -e "MGE_Type\tElement_Type\tWHOpriority">updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority1.txt
awk -F'\t' '{if($1!="MGE_Type"){str=$1"\t"$2;if(!(str in a)){a[str]=$8}else{a[str]=a[str]";"$8}}}END{for(str2 in a){if(a[str2]~/Critical priority/){print str2"\tCritical priority"}else if(a[str2]~/High priority/){print str2"\tHigh priority"}else if(a[str2]~/Medium priority/){print str2"\tMedium priority"}else{print str2"\t-"}}}' updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority0.txt>>updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority1.txt
cp updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority1.txt updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority.txt
awk -F'\t' 'NR==FNR{str=$1"\t"$2;a[str]=1}NR>FNR{str2=$1"\t"$2;if(!(str2 in a)){print str2"\tNone"}}' updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority1.txt updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02.txt|sort|uniq >>updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority.txt

rm updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority0.txt
rm updownstreamInfo/all_species_MGEresult_MobileElementFinder_DANMEL_BacAnt_allMGE_crossedLineage_genomeProp0.02_WHOpriority1.txt
#############################################################################################################################################################################################



############### 生成各风险评估体系各项指标结果 ################
cd ~/Feature_analysis
##沙门之前运行过mge的基因组的mge和argcombine的关联
sbatch -J run -A p_phage -p batch4 -n 10 -N 1 18_0_get_Salmon_prerunned_mgepath_run.sh
#4915102
awk -F'\t' 'NR==FNR{a[$1]=$0}NR>FNR{if($1 in a){print a[$1]}}' MGE_Salmonella_enterica_runedgenome_path.txt Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_nowhavemge.txt>MGE_Salmonella_enterica_runedgenome_path_2014use.txt
cat MGE_Salmonella_enterica_runedgenome_path_2014use.txt|while read i; do arr=($i); echo -e "perl $PWD/18_get_argAssomge_eachgenome.pl ~/Feature_analysis/ARGcombine_RefSeq_result/${arr[0]}/${arr[0]}_Salmonella_RefSeq_add_prokka_arg_fromdiamond_last_use_result.txt ${arr[1]}/${arr[0]}_MGE_merged_site.txt $PWD/MGE_prerunned/${arr[0]}_MGE_merged_ARGcombine_result.txt">>18_get_argAssomge_eachgenome_prerunnedmge_run.list;done
split -l 200 18_get_argAssomge_eachgenome_prerunnedmge_run.list rundata/18_get_argAssomge_eachgenome_prerunnedmge_run.list.split
ls rundata/18_get_argAssomge_eachgenome_prerunnedmge_run.list.split*|while read i
do
	echo "sh $PWD/${i}">>18_get_argAssomge_eachgenome_prerunnedmge_run.list.split.sh
done
sbatch -J mgemerge -A p_phage -p batch4 -n 10 -N 1 -a 1-373 ~/script/wrapper.sh 18_get_argAssomge_eachgenome_prerunnedmge_run.list.split.sh 1
#4915159

##沙门后补充mge的基因组的mge和argcombine的关联
##结果中没有DANMEL数据库比对的IS结果，前面融合时只融合了MobileElementFinder和BacAnt(发现14880自有的里面目前按照IS同家族关联的没有DANMEL的结果，公共病原项目之前运行的也没有DANMEL结果，因此本处没有再补充。如果后续要更换IS形成转座子的条件，再考虑把本次补充的10万数据加上DANMEL结果)
MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list;done
split -l 60 -d MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list rundata/MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list.split
ls rundata/MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list.split*|while read i; do echo -e "sh $PWD/$i">>MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list.split.sh;done
sed -n '1,1000p' MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list.split.sh>MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list.split_1.sh
sed -n '1001,$p' MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list.split.sh>MGE_salmone_bcrun_MobileElementFinder_DANMEL_BacAnt_run.list.split_2.sh


sed -n '2,$p' Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_sort_needrunmge.txt |cut -f1|while read i;do echo -e "perl $PWD/18_get_argAssomge_eachgenome.pl ~/Feature_analysis/ARGcombine_RefSeq_result/${i}/${i}_Salmonella_RefSeq_add_prokka_arg_fromdiamond_last_use_result.txt ~/Feature_analysis/MGE_bc_result/Salmonella_enterica/${i}/${i}_MGE_merged_site.txt ~/Feature_analysis/MGE_bc_result/Salmonella_enterica/${i}/${i}_MGE_merged_ARGcombine_result.txt">>18_get_argAssomge_eachgenome_bcmge_run.list;done
split -l 1000 18_get_argAssomge_eachgenome_bcmge_run.list rundata/18_get_argAssomge_eachgenome_bcmge_run.list.split
ls rundata/18_get_argAssomge_eachgenome_bcmge_run.list.split*|while read i
do
        echo "sh $PWD/${i}">>18_get_argAssomge_eachgenome_bcmge_run.list.split.sh
done
sbatch -J mgemerge -A p_phage -p batch2 -n 10 -N 1 -a 1-110 ~/script/wrapper.sh 18_get_argAssomge_eachgenome_bcmge_run.list.split.sh 1
#4930492
#4945980

#####合并沙门所用公共基因组的预警评分需要的信息
#ls use_topfeature|awk -F'_' '{print $1}'|grep -v "XGBoost"|grep -v "^pre"|while read i; do echo -e "sh $PWD/19_add_year_host_country_argAssomge.sh $i">>19_add_year_host_country_argAssomge_run.list;done
#sbatch -J merge -A p_phage -p batch3 -n 10 -N 1 -a 1-14 ~/script/wrapper.sh 19_add_year_host_country_argAssomge_run.list 1
##4916563
#目前只计算了已有mge的公共沙门基因组，还是10万的基因组待mge运行完后执行。


########### 自有14880基因组 mge关联 ####
sh 22_get_ownSalmon_argannoMatrix_run.list
cat ../Salmonella_MIC_usegenome_20240829.list |while read i; do echo -e "perl $PWD/23_get_own14880_argAssomge_eachgenome.pl ~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last.txt ~/MGE/MGE_result/Salmonella/${i}/${i}_MGE_merged_site.txt ~/MGE/MGE_result/Salmonella/${i}/${i}_MGE_merged_ARGcombine_result.txt">>23_get_own14880_argAssomge_eachgenome_run.list;done
split -l 100 23_get_own14880_argAssomge_eachgenome_run.list rundata/23_get_own14880_argAssomge_eachgenome_run.list.split
ls rundata/23_get_own14880_argAssomge_eachgenome_run.list.split*|while read i
do
        echo "sh $PWD/${i}">>23_get_own14880_argAssomge_eachgenome_run.list.split.sh
done
sbatch -J mgemerge -A p_phage -p batch2 -n 10 -N 1 -a 1-149 ~/script/wrapper.sh 23_get_own14880_argAssomge_eachgenome_run.list.split.sh 1
#4920716
#split -l 100 23_get_own14880_argAssomge_eachgenome_2_run.list rundata/23_get_own14880_argAssomge_eachgenome_2_run.list.split
#ls rundata/23_get_own14880_argAssomge_eachgenome_2_run.list.split*|while read i
#do
#        echo "sh $PWD/${i}">>23_get_own14880_argAssomge_eachgenome_2_run.list.split.sh
#done
#sbatch -J mgemerge -A p_phage -p batch2 -n 4 -N 1 -a 1-149 ~/script/wrapper.sh 23_get_own14880_argAssomge_eachgenome_2_run.list.split.sh 1
##4920865


##### 合并公共和自有每个基因组各注释到的特征的关联mge情况
sbatch -J mge -A p_phage -p batch4 -w node41 -n 20 -N 1 19_merge_184251_publice_mgeARGinfo_run.sh
#4930940
#4948335
sbatch -J mge -A p_phage -p batch4 -n 20 -N 1 19_merge_14880_own_mgeARGinfo_run.sh
#4930938
#4948358
## 上面特征仅考虑了基因，没考虑基因中SNV，因此补充上
sbatch -J mge -A p_phage -p batch3 -n 20 -N 1 19_merge_184251_publice_mgeARGinfo_SNV_run.sh
#5020093
#5126413
sbatch -J mge -A p_phage -p batch3 -n 20 -N 1 19_merge_14880_own_mgeARGinfo_SNV_run.sh
#5020094
#5126414

##### 合并公共和自有每个基因组各注释到的仅mge情况
sbatch -J mge -A p_phage -p batch4 -n 20 -N 1 20_merge_184251_publice_mgeinfo_run.sh
#4930942
#4948464
sbatch -J mge -A p_phage -p batch4 -n 20 -N 1 20_merge_14880_own_mgeinfo_run.sh
#4930944
#4948471
##### 合并公共和自有每个基因组各注释到的仅mge情况-添加mge在基因组位置信息
sbatch -J mge -A p_phage -p batch4 -w node40 -n 20 -N 1 20_merge_184251_publice_mgeinfo_addpos_run.sh
#4952492
sbatch -J mge -A p_phage -p batch4 -w node43 -n 20 -N 1 20_merge_14880_own_mgeinfo_addpos_run.sh
#4952493







####################################################################################################################################################################################
######################################  更换预警打分体系--去除沙门内部mge频率情况，改用arg+mge在输入的菌株基因组位置，改特征出现频率变化  ##########################################
####################################################################################################################################################################################
#### 生成每个基因组中特征+mge的位置（质粒、质粒+染色体、染色体，后期计算时质粒+染色体同质粒加分一致，均加1，染色体上不加分）
sbatch -J merge -A p_phage -p batch4 -n 20 -N 1 19_merge_14880_own_mgeARGinfo_addpos_run.sh 
#5288329
sbatch -J merge -A p_phage -p batch4 -n 20 -N 1 19_merge_14880_own_mgeARGinfo_SNV_addpos_run.sh
#5288331
sbatch -J merge -A p_phage -p batch2 -n 20 -N 1 19_merge_184251_publice_mgeARGinfo_addpos_run.sh
#5288332
sbatch -J merge -A p_phage -p batch4 -n 20 -N 1 19_merge_184251_publice_mgeARGinfo_SNV_addpos_run.sh
#5288333
#### combine ARG和SNV关联mge的结果
cp Salmonella_own_14880_mgeARGinfo_addpos.txt Salmonella_own_14880_mgeARGinfo_ARGaddSNV_addpos.txt
sed -n '2,$p' Salmonella_own_14880_mgeARGinfo_SNV_addpos.txt>>Salmonella_own_14880_mgeARGinfo_ARGaddSNV_addpos.txt
cp Salmonella_public_184251_mgeARGinfo_addpos.txt Salmonella_public_184251_mgeARGinfo_ARGaddSNV_addpos.txt
sed -n '2,$p' Salmonella_public_184251_mgeARGinfo_SNV_addpos.txt>>Salmonella_public_184251_mgeARGinfo_ARGaddSNV_addpos.txt
#### filt ANI & combine own and public
head -1 Salmonella_public_184251_mgeARGinfo_ARGaddSNV_addpos.txt>Salmonella_public_184251_mgeARGinfo_ARGaddSNV_addpos_filtANI96.txt
awk -F'\t' 'NR==FNR{a[$1]=1}NR>FNR{if($1 in a){print}}' Salmonella_public_184251_ANI96need_genomelist.txt Salmonella_public_184251_mgeARGinfo_ARGaddSNV_addpos.txt>>Salmonella_public_184251_mgeARGinfo_ARGaddSNV_addpos_filtANI96.txt
cp Salmonella_public_184251_mgeARGinfo_ARGaddSNV_addpos_filtANI96.txt Salmonella_own14880_public_184251_mgeARGinfo_ARGaddSNV_addpos_filtANI96.txt
sed -n '2,$p' Salmonella_own_14880_mgeARGinfo_ARGaddSNV_addpos.txt >>Salmonella_own14880_public_184251_mgeARGinfo_ARGaddSNV_addpos_filtANI96.txt

###### 生成每个基因组关联mge的结果
sbatch -J a -A p_phage -p batch4 -n 10 -N 1 -a 1-14 ~/script/wrapper.sh 29_2use_get_genome_usemge_info_run.list 1
#5288413
#有关联MGEs,但MGEs不跨种的CrossedLineage为”-“。

###### 筛选耐药的基因组
#perl 31_combine_predictRS.pl
ls predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_*_topfeature_matrix_add_mge.txt|sed 's#predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_##g'|sed 's#_topfeature_matrix_add_mge.txt##g'|while read i;do echo -e "perl $PWD/32_get_Rgenome.pl $i">>32_get_Rgenome_run.list;done
sbatch -J a -A p_phage -p batch4 -n 10 -N 1 -a 1-14 ~/script/wrapper.sh 32_get_Rgenome_run.list 1
#5288447

########## 每种抗生素各耐药基因组风险打分
### 特征趋势打分改为最低为正值
mv predict_use_data/Antibiotic_features_CAGRscore_2017-2022_thread5.txt predict_use_data/Antibiotic_features_CAGRscore_2017-2022_thread5_pre.txt
awk -F'\t' '{if($1=="Feature"){print }else{if($2=="0"){print $1"\t1"}else if($2~/^-/){print $1"\t0"}else{print $0}}}' predict_use_data/Antibiotic_features_CAGRscore_2017-2022_thread5_pre.txt>predict_use_data/Antibiotic_features_CAGRscore_2017-2022_thread5.txt

sbatch -J risk -A p_phage -p batch4 -n 10 -N 1 -a 1-14 ~/script/wrapper.sh 33_risk_score_run.list 1
#5288538
#5291457
#### 整合各类抗生素打分
perl 33_2_combine_risk_score.pl
### 打分取整
awk -F'\t' '{if($1=="Genome"){print $0}else{printf $1;for(i=2;i<=NF;i++){if($i=="-"){printf "\t-"}else{printf "\t%.f",$i}}printf "\n"}}' predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR.txt>predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_int.txt

### 整合时间国家和打分结果
awk -F'\t' 'NR==FNR{a[$1]=$(NF-2)"\t"$NF}NR>FNR{if($1 in a){print a[$1]"\t"$0}}' predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_GEN_topfeature_matrix_add_year_host_country.txt predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR.txt>predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_add_year_country.txt
### 划分各个基因组的风险等级
awk -F'\t' '{if($1=="Genome"){printf $0}else{printf "\n"$1;for(i=2;i<=NF;i++){if($i=="-"){printf "\t-"}else if($i>=12){printf "\tI"}else if($i>=4){printf "\tII"}else{printf "\tIII"}}}}' predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR.txt>predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel.txt
awk -F'\t' 'NR==FNR{a[$1]=$(NF-2)"\t"$NF}NR>FNR{if($1 in a){print a[$1]"\t"$0}}' predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_GEN_topfeature_matrix_add_year_host_country.txt predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel.txt>predict_use_data/All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel_add_year_country.txt



#### ownuse14846_publicANI96_173348的基因组预测表型结果，添加分离年份和宿主、国家信息
awk -F'\t' 'NR==FNR{if($1!="Genome"){a[$1]=$(NF-2)"\t"$(NF-1)"\t"$NF}}NR>FNR{if($1=="Genome"){print $0"\tyear\thost\tcountry"}else{if($1 in a){print $0"\t"a[$1]}}}' predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_AMP_topfeature_matrix_add_year_host_country.txt predict_use_data/Salmonella_enterica_ownuse14880_public_184252_Allantibiotic_topfeature_matrix_result.txt>predict_use_data/Salmonella_enterica_ownuse14846_publicANI96_173348_Allantibiotic_topfeature_matrix_result_add_year_host_country.txt

