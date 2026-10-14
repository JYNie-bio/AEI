#!/bin/sh
cd ~/ARG_MIC
antibiotic="CHL"
awk -F'\t' -v antibiotic=$antibiotic 'NR==FNR{if($1==antibiotic && $NF=="near") a[$2]=1;next}FNR==1{print; next}($1 in a){print}' ../chewbbaca/test_near_vs_independent_${antibiotic}.txt RS_811/14880file_ARGanno_Protein_result_copynumber_addPhenotype_${antibiotic}_verify_test.txt>RS_811/14880file_ARGanno_Protein_result_copynumber_addPhenotype_${antibiotic}_verify_test_cgmlst_near_inOwnValid.txt
awk -F'\t' -v antibiotic=$antibiotic 'NR==FNR{if($1==antibiotic && $NF=="independent") a[$2]=1;next}FNR==1{print; next}($1 in a){print}' ../chewbbaca/test_near_vs_independent_${antibiotic}.txt RS_811/14880file_ARGanno_Protein_result_copynumber_addPhenotype_${antibiotic}_verify_test.txt>RS_811/14880file_ARGanno_Protein_result_copynumber_addPhenotype_${antibiotic}_verify_test_cgmlst_independent_inOwnValid.txt
antibiotic="AMP"
antibiotic="AMS"
antibiotic="CTX"
antibiotic="TET"
awk -F'\t' -v antibiotic=$antibiotic 'NR==FNR{if($1==antibiotic && $NF=="near") a[$2]=1;next}FNR==1{print; next}($1 in a){print}' ../chewbbaca/test_near_vs_independent_${antibiotic}.txt RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_${antibiotic}_verify_test.txt>RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_${antibiotic}_verify_test_cgmlst_near_inOwnValid.txt
awk -F'\t' -v antibiotic=$antibiotic 'NR==FNR{if($1==antibiotic && $NF=="independent") a[$2]=1;next}FNR==1{print; next}($1 in a){print}' ../chewbbaca/test_near_vs_independent_${antibiotic}.txt RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_${antibiotic}_verify_test.txt>RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_${antibiotic}_verify_test_cgmlst_independent_inOwnValid.txt

