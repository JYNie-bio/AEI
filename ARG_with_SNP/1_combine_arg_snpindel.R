args=commandArgs(T)

setwd("~/ARG_with_SNP")
###### ARG use copynumber
argdata<-read.table("~/RefSeq_Anno/14880file_ARGanno_Protein_result_copynumber_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata1<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata3<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_YesNoNA_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)

snpdata1_1<-snpdata1[rownames(argdata),18:ncol(snpdata1)]
snpdata1_2<-ifelse(snpdata1_1 == 1, "Yes", "No")
snpdata3_2<-snpdata3[rownames(argdata),18:ncol(snpdata3)]

alldata1<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata1_2)
alldata3<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata3_2)
colnames(alldata1)[1]<-"Genome_name"
colnames(alldata3)[1]<-"Genome_name"

write.table(alldata1,"14880file_ARGanno_copynumber_SNPmatrix959090_absence01_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)
write.table(alldata3,"14880file_ARGanno_copynumber_SNPmatrix959090_YesNoNA_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)

###### ARG use absence01
argdata<-read.table("~/RefSeq_Anno/14880file_ARGanno_Protein_result_absence01_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata1<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata3<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_YesNoNA_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)

snpdata1_2<-snpdata1[rownames(argdata),18:ncol(snpdata1)]
snpdata3_2<-snpdata3[rownames(argdata),18:ncol(snpdata3)]

alldata1<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata1_2)
alldata3<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata3_2)
colnames(alldata1)[1]<-"Genome_name"
colnames(alldata3)[1]<-"Genome_name"

write.table(alldata1,"14880file_ARGanno_absence01_SNPmatrix959090_absence01_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)
write.table(alldata3,"14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)

###### ARG use copynumber
argdata<-read.table("~/RefSeq_Anno/14880file_ARGanno_Protein_result_copynumber_add_mutation_models_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata1<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_add_mutation_models_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata3<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_YesNoNA_add_mutation_models_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)

snpdata1_1<-snpdata1[rownames(argdata),18:ncol(snpdata1)]
snpdata1_2<-ifelse(snpdata1_1 == 1, "Yes", "No")
snpdata3_2<-snpdata3[rownames(argdata),18:ncol(snpdata3)]

alldata1<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata1_2)
alldata3<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata3_2)
colnames(alldata1)[1]<-"Genome_name"
colnames(alldata3)[1]<-"Genome_name"

write.table(alldata1,"14880file_ARGanno_copynumber_SNPmatrix959090_absence01_add_mutation_models_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)
write.table(alldata3,"14880file_ARGanno_copynumber_SNPmatrix959090_YesNoNA_add_mutation_models_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)

###### ARG use absence01
argdata<-read.table("~/RefSeq_Anno/14880file_ARGanno_Protein_result_absence01_add_mutation_models_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata1<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_add_mutation_models_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)
snpdata3<-read.table("~/RefSeq_Anno/RefSeq_Anno_Result_14880file_addinfo_add_prokka_arg_fromdiamond_havearg_last_use_result_last_SNPmatrix959090_YesNoNA_add_mutation_models_addPhenotype.txt",header=T,quote='',sep="\t",check.names = F,row.names=1)

snpdata1_2<-snpdata1[rownames(argdata),18:ncol(snpdata1)]
snpdata3_2<-snpdata3[rownames(argdata),18:ncol(snpdata3)]

alldata1<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata1_2)
alldata3<-cbind(rownames(argdata),argdata[,1:(ncol(argdata)-1)],snpdata3_2)
colnames(alldata1)[1]<-"Genome_name"
colnames(alldata3)[1]<-"Genome_name"

write.table(alldata1,"14880file_ARGanno_absence01_SNPmatrix959090_absence01_add_mutation_models_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)
write.table(alldata3,"14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_add_mutation_models_addPhenotype.txt",quote=F,sep="\t",col.names=T,row.names=F)

