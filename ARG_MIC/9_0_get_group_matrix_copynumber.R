args=commandArgs(T)
#args is AMP
ptm <- proc.time()
print(ptm)


setwd("~/ARG_MIC/RS_811")

#### train\val\test
datatest<-read.table(paste0("~/chewbbaca/split_results/",args[1],"_test.txt"))
data_all<-read.table("~/RefSeq_Anno/14880file_ARGanno_Protein_result_copynumber_addPhenotype.txt",header = T,sep="\t",quote = '',comment.char = '',check.names = F,row.names=1,stringsAsFactors = F)
data_arg<-data_all[,c(17:ncol(data_all),which(colnames(data_all)==args[1]))]
colnames(data_arg)[ncol(data_arg)]<-"Group"

datatest<-read.table(paste0("~/chewbbaca/split_results/",args[1],"_test.txt"),header = T,sep="\t",quote = '',comment.char = '',check.names = F,stringsAsFactors = F)
data_arg_test<-data_arg[datatest[,1],]
write.table(data_arg_test,paste0("14880file_ARGanno_Protein_result_copynumber_addPhenotype_",args[1],"_verify_test_stcgmlst.txt"),sep="\t",quote=F,row.names = T,col.names = T)


for(i in 1:9){
	datatrain<-read.table(paste0("~/chewbbaca/split_results/",args[1],"_fold0",i,"_train.txt"),header = T,sep="\t",quote = '',comment.char = '',check.names = F,stringsAsFactors = F)
	dataval<-read.table(paste0("~/chewbbaca/split_results/",args[1],"_fold0",i,"_val.txt"),header = T,sep="\t",quote = '',comment.char = '',check.names = F,stringsAsFactors = F)
	data_arg_train<-data_arg[datatrain[,1],]
	data_arg_val<-data_arg[dataval[,1],]
	data_arg_train$Group<-gsub("R","1",data_arg_train$Group)
	data_arg_train$Group<-gsub("S","0",data_arg_train$Group)
	data_arg_val$Group<-gsub("R","1",data_arg_val$Group)
        data_arg_val$Group<-gsub("S","0",data_arg_val$Group)
	write.table(data_arg_train,paste0("~/ARG_MIC/AllFeature/XGBoost_ARG_copynumber_fold_",i,"_train_samples_10fold_",args[1],"_stcgmlst.txt"),sep="\t",quote=F,row.names = T,col.names = T)
	write.table(data_arg_val,paste0("~/ARG_MIC/AllFeature/XGBoost_ARG_copynumber_fold_",i,"_val_samples_10fold_",args[1],"_stcgmlst.txt"),sep="\t",quote=F,row.names = T,col.names = T)
}

datatrain<-read.table(paste0("~/chewbbaca/split_results/",args[1],"_fold10_train.txt"),header = T,sep="\t",quote = '',comment.char = '',check.names = F,stringsAsFactors = F)
dataval<-read.table(paste0("~/chewbbaca/split_results/",args[1],"_fold10_val.txt"),header = T,sep="\t",quote ='',comment.char = '',check.names = F,stringsAsFactors = F)
datatrainalllist<-c(datatrain[,1],dataval[,1])
data_arg_trainall<-data_arg[datatrainalllist,]
write.table(data_arg_trainall,paste0("14880file_ARGanno_Protein_result_copynumber_addPhenotype_",args[1],"_train_stcgmlst.txt"),sep="\t",quote=F,row.names = T,col.names = T)
data_arg_train<-data_arg[datatrain[,1],]
data_arg_val<-data_arg[dataval[,1],]
data_arg_train$Group<-gsub("R","1",data_arg_train$Group)
data_arg_train$Group<-gsub("S","0",data_arg_train$Group)
data_arg_val$Group<-gsub("R","1",data_arg_val$Group)
data_arg_val$Group<-gsub("S","0",data_arg_val$Group)
write.table(data_arg_train,paste0("~/ARG_MIC/AllFeature/XGBoost_ARG_copynumber_fold_10_train_samples_10fold_",args[1],"_stcgmlst.txt"),sep="\t",quote=F,row.names = T,col.names = T)
write.table(data_arg_val,paste0("~/ARG_MIC/AllFeature/XGBoost_ARG_copynumber_fold_10_val_samples_10fold_",args[1],"_stcgmlst.txt"),sep="\t",quote=F,row.names = T,col.names = T)


ptm2<-proc.time() - ptm
print(ptm2)
