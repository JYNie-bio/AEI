#!/bin/bash
#Mobile genetic elements identifying
#Description:This program can automatic identify MGEs and merge them into two files.
#use as: sh MGE_gcPathogen_bacteria_pipeline.sh -O ~/MGE_result/ -t 20 -I /genome/GCA_002787645.1.fna -i 90 -s 80 -q 80 -e 0.00001
#use as: If do not provide gbk and gff3 file, this pipeline will automatic generate.
#version="1.0 version";
#echo -e "*************\n*$version*\n*************";

#############################################################################################################################################################
#NOTE: the contig name must begain with a letter, not a number. Moreover, contig name cannot contain a colon, otherwise, the result will be affected.
#############################################################################################################################################################

while getopts ":O:I:t:g:f:k:S:i:q:s:e:" opt ;do
        case $opt in
		O) outdir=$OPTARG;; ###output directory
                I) inputfafile=$OPTARG;; ###input fasta or fna or fa file, including the absolute path to the location.
                t) threads=$OPTARG;; ###run threads
		g) gbkfile=$OPTARG;; ###genbank format file for input fasta file, including the absolute path to the location.
		f) gfffile=$OPTARG;; ###gff3 format file for input fasta file, including the absolute path to the location.
		k) fakingdom=$OPTARG;; ###Kingdom to which the fasta file belongs.
		S) species0=$OPTARG;; ###Species for genome
		i) identity=$OPTARG;; ###identity
		q) qcov=$OPTARG;; ###query coverage
		s) scov=$OPTARG;; ###subject coverage
		e) evalue=$OPTARG;; ###evalue
        esac
done

if [[ ! $outdir =~ .*\/$ ]]
then
        outdir=$outdir"/"
fi

identity_c=`echo "scale=1;$identity/100"|bc`
qcov_c=`echo "scale=1;$qcov/100"|bc`
scov_c=`echo "scale=1;$scov/100"|bc`

if [[ $inputfafile =~ .*gz$ ]]
then
        if [[ $inputfafile =~ .*fna.gz$ ]]
        then
                filearr=(${inputfafile//\// })
                filename=${filearr[${#filearr[*]}-1]}
                species=${filearr[${#filearr[*]}-2]}
                faprefix=${filename%.fna.gz}
                fasuffix="fna"
                mkdir -p ${outdir}${faprefix}/temp
                cp ${inputfafile} ${outdir}${faprefix}/temp
                gunzip ${outdir}${faprefix}/temp/${filename}
                inputfafile=${outdir}${faprefix}/temp/${faprefix}.${fasuffix}
		cd ${outdir}${faprefix}/temp
		sed -i 's#\r##g' ${inputfafile}
                sed -i "s#:#_#g" ${inputfafile}
		. ~/miniconda3/bin/activate samtools1.21
		samtools faidx ${inputfafile}
		conda deactivate
        elif [[ $inputfafile =~ .*fasta.gz$ ]]
        then
                filearr=(${inputfafile//\// })
                filename=${filearr[${#filearr[*]}-1]}
                species=${filearr[${#filearr[*]}-2]}
                faprefix=${filename%.fasta.gz}
                fasuffix="fasta"
                mkdir -p ${outdir}${faprefix}/temp
                cp ${inputfafile} ${outdir}${faprefix}/temp
                gunzip ${outdir}${faprefix}/temp/${filename}
                inputfafile=${outdir}${faprefix}/temp/${faprefix}.${fasuffix}
                cd ${outdir}${faprefix}/temp
		sed -i 's#\r##g' ${inputfafile}
                sed -i "s#:#_#g" ${inputfafile}
                . ~/miniconda3/bin/activate samtools1.21
                samtools faidx ${inputfafile}
                conda deactivate
        elif [[ $inputfafile =~ .*fa.gz$ ]]
        then
                filearr=(${inputfafile//\// })
                filename=${filearr[${#filearr[*]}-1]}
                species=${filearr[${#filearr[*]}-2]}
                faprefix=${filename%.fa.gz}
                fasuffix="fa"
                mkdir -p ${outdir}${faprefix}/temp
                cp ${inputfafile} ${outdir}${faprefix}/temp
                gunzip ${outdir}${faprefix}/temp/${filename}
                inputfafile=${outdir}${faprefix}/temp/${faprefix}.${fasuffix}
                cd ${outdir}${faprefix}/temp
		sed -i 's#\r##g' ${inputfafile}
                sed -i "s#:#_#g" ${inputfafile}
                . ~/miniconda3/bin/activate samtools1.21
                samtools faidx ${inputfafile}
                conda deactivate
        else
                echo "ERROR: The file is not fasta.gz or fna.gz or fa.gz ending file.">>${outdir}${faprefix}/pipeline.log
                echo "ERROR: The file is not fasta.gz or fna.gz or fa.gz ending file."
        fi
	isescanspecies="temp"
elif [[ $inputfafile =~ .*fna$ ]]
then
        filearr=(${inputfafile//\// })
        filename=${filearr[${#filearr[*]}-1]}
	species=${filearr[${#filearr[*]}-2]}
        faprefix=${filename%.fna}
        fasuffix="fna"
	mkdir -p ${outdir}${faprefix}/temp
	cp ${inputfafile} ${outdir}${faprefix}/temp
        cd ${outdir}${faprefix}/temp
	sed -i 's#\r##g' ${filename}
        sed -i "s#:#_#g" ${filename}
        . ~/miniconda3/bin/activate samtools1.21
        samtools faidx ${filename}
        conda deactivate
	inputfafile=${outdir}${faprefix}/temp/${faprefix}.${fasuffix}
	isescanspecies="temp"
elif [[ $inputfafile =~ .*fasta$ ]]
then
        filearr=(${inputfafile//\// })
        filename=${filearr[${#filearr[*]}-1]}
	species=${filearr[${#filearr[*]}-2]}
        faprefix=${filename%.fasta}
        fasuffix="fasta"
        mkdir -p ${outdir}${faprefix}/temp
        cp ${inputfafile} ${outdir}${faprefix}/temp
        cd ${outdir}${faprefix}/temp
	sed -i 's#\r##g' ${filename}
        sed -i "s#:#_#g" ${filename}
        . ~/miniconda3/bin/activate samtools1.21
        samtools faidx ${filename}
        conda deactivate
	inputfafile=${outdir}${faprefix}/temp/${faprefix}.${fasuffix}
	isescanspecies="temp"
elif [[ $inputfafile =~ .*fa$ ]]
then
        filearr=(${inputfafile//\// })
        filename=${filearr[${#filearr[*]}-1]}
	species=${filearr[${#filearr[*]}-2]}
        faprefix=${filename%.fa}
        fasuffix="fa"
        mkdir -p ${outdir}${faprefix}/temp
        cp ${inputfafile} ${outdir}${faprefix}/temp
        cd ${outdir}${faprefix}/temp
	sed -i 's#\r##g' ${filename}
        sed -i "s#:#_#g" ${filename}
        . ~/miniconda3/bin/activate samtools1.21
        samtools faidx ${filename}
        conda deactivate
	inputfafile=${outdir}${faprefix}/temp/${faprefix}.${fasuffix}
	isescanspecies="temp"
else
        echo "ERROR: The file is not fasta or fna or fa ending file.">>${outdir}${faprefix}/pipeline.log
        echo "ERROR: The file is not fasta or fna or fa ending file."
	exit 2
fi

if [ -n "$species0" ];then
	species=${species0}
fi

if [ ! -n "$outdir" ] || [ ! -n "$inputfafile" ] || [ ! -n "$identity" ] || [ ! -n "$qcov" ] || [ ! -n "$scov" ] || [ ! -n "$evalue" ];then
	echo "************************************************************************
	ERROR: Missing required parameter.
	Usage: sh MGE_gcPathogen_bacteria_pipeline.sh -O outdirectory -I inputfastafile_addpath -i identity -q qcov -s scov -e evalue
************************************************************************">>${outdir}${faprefix}/pipeline.log
	echo "************************************************************************
        ERROR: Missing required parameter.
        Usage: sh MGE_gcPathogen_bacteria_pipeline.sh -O outdirectory -I inputfastafile_addpath -i identity -q qcov -s scov -e evalue
************************************************************************"
	exit 2
fi

#gfffilearr=(${gfffile//\// })
#gffspecies=${gfffilearr[${#gfffilearr[*]}-2]}

if [ -n "$fakingdom" ];then
	if [ "$fakingdom" != "Archaea" ] && [ "$fakingdom" != "Bacteria" ]
	then
		echo "ERROR: The fasta kingdom is wrong. Need to check fakingdom and rerun this program.">>${outdir}${faprefix}/pipeline.log
		echo "ERROR: The fasta kingdom is wrong. Need to check fakingdom and rerun this program."
		exit 2
	fi
else
	fakingdom="Bacteria"
fi

if [ ! -n "$threads" ];then
	if [[ -n "$SLURM_NTASKS" ]];then
	        threads=$SLURM_NTASKS
	else
		threads=10
	fi
fi


#### identity MGE ####
if [ ! -d "$outdir" ];then
	mkdir -p $outdir
fi
mkdir -p ${outdir}${faprefix}/plasmidfinder
mkdir -p ${outdir}${faprefix}/genomad
mkdir -p ${outdir}${faprefix}/MobileElementFinder


cd $outdir
OLDPATH=$PATH

date 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log


ProkkaStartTime=`date +'%Y-%m-%d %H:%M:%S'`
ProkkaStartTime_s=$(date --date="$ProkkaStartTime" +%s)
##Run prokka and rewrite gff file
if [ ! -n "$gbkfile" ] && [ ! -n "$gfffile" ];then
	echo "############# prokka ####################"
	echo "############# prokka ####################">>${outdir}${faprefix}/pipeline.log
        . ~/miniconda3/bin/activate prokka
        prokka --kingdom ${fakingdom} --outdir ${outdir}${faprefix}/prokka --prefix ${faprefix} --cpus ${threads} ${inputfafile} --gcode 11 --force 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
        conda deactivate
	gbkfile=${outdir}${faprefix}/prokka/${faprefix}.gbk
	gfffile=${outdir}${faprefix}/prokka/${faprefix}.gff
	echo "We will use prokka result for gbk and gff"
        echo "We will use prokka result for gbk and gff">>${outdir}${faprefix}/pipeline.log
elif [ -n "$gbkfile" ] && [ ! -n "$gfffile" ];then
	echo "############# prokka ####################"
	echo "############# prokka ####################">>${outdir}${faprefix}/pipeline.log
        . ~/miniconda3/bin/activate prokka
        prokka --kingdom ${fakingdom} --outdir ${outdir}${faprefix}/prokka --prefix ${faprefix} --cpus ${threads} ${inputfafile} --gcode 11 --force 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
        conda deactivate
        gfffile=${outdir}${faprefix}/prokka/${faprefix}.gff
	echo "We will use prokka result for gff"
	echo "We will use prokka result for gff">>${outdir}${faprefix}/pipeline.log
elif [ ! -n "$gbkfile" ] && [ -n "$gfffile" ];then
	echo "############# prokka ####################"
	echo "############# prokka ####################">>${outdir}${faprefix}/pipeline.log
        . ~/miniconda3/bin/activate prokka
        prokka --kingdom ${fakingdom} --outdir ${outdir}${faprefix}/prokka --prefix ${faprefix} --cpus ${threads} ${inputfafile} --gcode 11 --force 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
        conda deactivate
        gbkfile=${outdir}${faprefix}/prokka/${faprefix}.gbk
	echo "We will use prokka result for gbk"
        echo "We will use prokka result for gbk">>${outdir}${faprefix}/pipeline.log
fi

perl ~/script/prokkagff_repair.pl ${gfffile} ${gbkfile} ${gbkfile}2
mv ${gbkfile}2 ${gbkfile}
sed -i 's#*$##g' ${outdir}${faprefix}/prokka/${faprefix}.faa

ProkkaEndTime=`date +'%Y-%m-%d %H:%M:%S'`
ProkkaEndTime_s=$(date --date="$ProkkaEndTime" +%s)
echo "Prokka run time：$ProkkaStartTime ---> $ProkkaEndTime "$((ProkkaEndTime_s-ProkkaStartTime_s))"s" >>${outdir}${faprefix}/pipeline.log



########################################################################################################
##plasmid-plasmidfinder
plasmidfinderStartTime=`date +'%Y-%m-%d %H:%M:%S'`
plasmidfinderStartTime_s=$(date --date="$plasmidfinderStartTime" +%s)
date 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
. ~/miniconda3/bin/activate plasmidfinder
echo "############# plasmidfinder ####################"
echo "############# plasmidfinder ####################">>${outdir}${faprefix}/pipeline.log
plasmidfinder.py -i ${inputfafile} -o ${outdir}${faprefix}/plasmidfinder --extented_output 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
temp0=`grep "KeyError" ${outdir}${faprefix}/pipeline.log`
if [ ! -z "$temp0" ];then
        sed -i 's#KeyError#plasmidfinderERROR#g' ${outdir}${faprefix}/pipeline.log
fi
conda deactivate
plasmidfinderEndTime=`date +'%Y-%m-%d %H:%M:%S'`
plasmidfinderEndTime_s=$(date --date="$plasmidfinderEndTime" +%s)
echo "plasmidfinder run time：$plasmidfinderStartTime ---> $plasmidfinderEndTime "$((plasmidfinderEndTime_s-plasmidfinderStartTime_s))"s" >>${outdir}${faprefix}/pipeline.log

##plasmid-platon
PlatonStartTime=`date +'%Y-%m-%d %H:%M:%S'`
PlatonStartTime_s=$(date --date="$PlatonStartTime" +%s)
date 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
. ~/miniconda3/bin/activate platon1.7
echo "############# platon ####################"
echo "############# platon ####################">>${outdir}${faprefix}/pipeline.log
platon --db ~/software_test/platon/db --output ${outdir}${faprefix}/platon --verbose --threads ${threads} ${inputfafile}
conda deactivate
PlatonEndTime=`date +'%Y-%m-%d %H:%M:%S'`
PlatonEndTime_s=$(date --date="$PlatonEndTime" +%s)
echo "Platon run time：$PlatonStartTime ---> $PlatonEndTime "$((PlatonEndTime_s-PlatonStartTime_s))"s" >>${outdir}${faprefix}/pipeline.log

##plasmid-plasmer
PlasmerStartTime=`date +'%Y-%m-%d %H:%M:%S'`
PlasmerStartTime_s=$(date --date="$PlatonStartTime" +%s)
date 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
echo "############# plasmer ####################"
echo "############# plasmer ####################">>${outdir}${faprefix}/pipeline.log
export PATH=~/miniconda3/envs/plasmer/bin:$PATH
Plasmer -g  ${inputfafile} -p ${faprefix} -d ~/database/plasmer_db -t ${threads} -o ${outdir}${faprefix}/plasmer
PlasmerEndTime=`date +'%Y-%m-%d %H:%M:%S'`
PlasmerEndTime_s=$(date --date="$PlasmerEndTime" +%s)
echo "Plasmer run time：$PlasmerStartTime ---> $PlasmerEndTime "$((PlasmerEndTime_s-PlasmerStartTime_s))"s" >>${outdir}${faprefix}/pipeline.log

##genomad
geNomadStartTime=`date +'%Y-%m-%d %H:%M:%S'`
geNomadStartTime_s=$(date --date="$PlatonStartTime" +%s)
date 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
echo "############# geNomad ####################"
echo "############# geNomad ####################">>${outdir}${faprefix}/pipeline.log
. ~/miniconda3/bin/activate genomad
genomad end-to-end --cleanup -t $threads ${inputfafile} ${outdir}${faprefix}/genomad ~/databases/geNomad/genomad_db
geNomadEndTime=`date +'%Y-%m-%d %H:%M:%S'`
geNomadEndTime_s=$(date --date="$geNomadEndTime" +%s)
echo "geNomad run time：$geNomadStartTime ---> $geNomadEndTime "$((geNomadEndTime_s-geNomadStartTime_s))"s" >>${outdir}${faprefix}/pipeline.log


##transposon & integron
BacAntStartTime=`date +'%Y-%m-%d %H:%M:%S'`
BacAntStartTime_s=$(date --date="$BacAntStartTime" +%s)
date 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
. ~/miniconda3/bin/activate bacant
echo "############# bacant ####################"
echo "############# bacant ####################">>${outdir}${faprefix}/pipeline.log
export PATH=~/software/ncbi-blast-2.7.1+/bin:$PATH
bacant -n ${inputfafile} -c $scov,$scov,$scov -i $identity,$identity,$identity -o ${outdir}${faprefix}/bacant -t ${threads} 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
conda deactivate
PATH=$OLDPATH
BacAntEndTime=`date +'%Y-%m-%d %H:%M:%S'`
BacAntEndTime_s=$(date --date="$BacAntEndTime" +%s)
echo "BacAnt run time：$BacAntStartTime ---> $BacAntEndTime "$((BacAntEndTime_s-BacAntStartTime_s))"s" >>${outdir}${faprefix}/pipeline.log

##MobileElementFinder
MobileElementFinderStartTime=`date +'%Y-%m-%d %H:%M:%S'`
MobileElementFinderStartTime_s=$(date --date="$MobileElementFinderStartTime" +%s)
cd $outdir
echo "############# MobileElementFinder ####################"
echo "############# MobileElementFinder ####################">>${outdir}${faprefix}/pipeline.log
. ~/miniconda3/bin/activate MobileElementFinder1.1.2
mkdir -p ${outdir}${faprefix}/temp_Mobile
mefinder find --temp-dir ${outdir}${faprefix}/temp_Mobile --contig ${outdir}${faprefix}/temp/${faprefix}.${fasuffix} -t ${threads} --min-coverage $scov_c ${outdir}${faprefix}/MobileElementFinder/${faprefix}.MobileElementFinder.out 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
rm -r ${outdir}${faprefix}/temp_Mobile
tempa=`grep "subprocess.CalledProcessError: Command" ${outdir}${faprefix}/pipeline.log`
if [ ! -z "$tempa" ];then
        sed -i 's#CalledProcessError##g' ${outdir}${faprefix}/pipeline.log
        sed -i 's#ExitCodeError##g' ${outdir}${faprefix}/pipeline.log
        sed -i 's#BlastDatabaseError##g' ${outdir}${faprefix}/pipeline.log
        sed -i 's#me_finder.errors##g' ${outdir}${faprefix}/pipeline.log
        perl ~/script/fastafile_processMultiATCG.pl ${outdir}${faprefix}/temp/${faprefix}.${fasuffix} ${outdir}${faprefix}/temp/${faprefix}.${fasuffix}2
	mkdir -p ${outdir}${faprefix}/temp_Mobile
        mefinder find --temp-dir ${outdir}${faprefix}/temp_Mobile --contig ${outdir}${faprefix}/temp/${faprefix}.${fasuffix}2 -t ${threads} --min-coverage $scov_c ${outdir}${faprefix}/MobileElementFinder/${faprefix}.MobileElementFinder.out 2>&1 |tee -a ${outdir}${faprefix}/pipeline.log
	rm -r ${outdir}${faprefix}/temp_Mobile
fi
tempb=`grep "KeyError" ${outdir}${faprefix}/pipeline.log`
if [ ! -z "$tempb" ];then
	sed -i 's#KeyError##g' ${outdir}${faprefix}/pipeline.log
	perl ~/script/fastafile_processMultiATCG.pl ${outdir}${faprefix}/temp/${faprefix}.${fasuffix} ${outdir}${faprefix}/temp/${faprefix}.${fasuffix}2
	mkdir -p ${outdir}${faprefix}/temp_Mobile
	mefinder find --temp-dir ${outdir}${faprefix}/temp_Mobile --contig ${outdir}${faprefix}/temp/${faprefix}.${fasuffix}2 -t ${threads} --min-coverage $scov_c ${outdir}${faprefix}/MobileElementFinder/${faprefix}.MobileElementFinder.out 2>&1 |tee -a ${outdir}${faprefix}/pipeline.log
	rm -r ${outdir}${faprefix}/temp_Mobile
fi
conda deactivate
MobileElementFinderEndTime=`date +'%Y-%m-%d %H:%M:%S'`
MobileElementFinderEndTime_s=$(date --date="$MobileElementFinderEndTime" +%s)
echo "MobileElementFinder run time：$MobileElementFinderStartTime ---> $MobileElementFinderEndTime "$((MobileElementFinderEndTime_s-MobileElementFinderStartTime_s))"s" >>${outdir}${faprefix}/pipeline.log



#######################################################################################################
##### deal result #####
## deal with platon result
echo "############# deal platon result ####################">>${outdir}${faprefix}/pipeline.log
echo "############# deal platon result ####################"
perl ~/script/platon_result_proceed2.pl ~/database/ncbi_plasmid_download/plasmids_ncbi_20241010.txt ${outdir}${faprefix}/platon/${faprefix}.tsv ${outdir}${faprefix}/platon/${faprefix}.json ${outdir}${faprefix}/temp/platon_result_proceed.txt
## deal with plasmidfinder result
echo "############# deal plasmidfinder result ####################">>${outdir}${faprefix}/pipeline.log
echo "############# deal plasmidfinder result ####################"
perl ~/script/plasmidfinder_process_custom_cov_identity.pl ${outdir} ${faprefix} ${fasuffix} ${outdir}${faprefix}/temp/ plasmidfinder_result_proceed_custom_cov_identity.txt $identity $scov $qcov 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
perl ~/script/plasmidfinder_process_custom_cov_identity.pl ${outdir} ${faprefix} ${fasuffix} ${outdir}${faprefix}/temp/ plasmidfinder_result_proceed.txt 90 80 80 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log


####################################################################################################################################
#### merge all MGE ####
### result file is in each sample folder, named MGE_merged_site.txt and MGE_merged_num.txt ####
echo "############# MGEmerge ####################">>${outdir}${faprefix}/pipeline.log
echo "############# MGEmerge ####################"
perl ~/script/MGE_merge14_custom_cov_identity.pl ${outdir} ${faprefix} ${isescanspecies} ${fasuffix} ${faprefix}_MGE_merged_site_custom_cov_identity.txt $qcov_c $scov_c $identity_c $evalue 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
### filt 908080
perl ~/script/MGE_merge14_908080.pl ${outdir} ${faprefix} ${isescanspecies} ${fasuffix} ${faprefix}_MGE_merged_site.txt 0.8 0.8 0.9 0.00001 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log

echo "############# MGEmerge_allcds ####################">>${outdir}${faprefix}/pipeline.log
echo "########## MGEmerge_allcds ###############";
perl ~/script/get_MGEmerged_eachsoftware_prokkacds14.pl ${gfffile} ${outdir}${faprefix}/ ${faprefix}_MGE_merged_site_custom_cov_identity.txt ${faprefix}_MGE_merged_cds_result_custom_cov_identity.txt 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log
### 908080MGEresult merge allcds
perl ~/script/get_MGEmerged_eachsoftware_prokkacds14.pl ${gfffile} ${outdir}${faprefix}/ ${faprefix}_MGE_merged_site.txt ${faprefix}_MGE_merged_cds_result.txt 2>&1 | tee -a ${outdir}${faprefix}/pipeline.log


rm ${outdir}${faprefix}/temp/${faprefix}.${fasuffix}
rm ${outdir}${faprefix}/temp/${faprefix}.${fasuffix}2

