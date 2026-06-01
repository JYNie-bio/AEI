args=commandArgs(T)
num=as.numeric(args[1])
#num is 2:16
ptm <- proc.time()
print(ptm)


library(Boruta)
setwd("~/ARG_with_SNP/RS_811")
tempstr="";
if(args[2]=="1"){
  tempstr="_add_mutation_models"
}

#### 8:2 split train\test
data_all<-read.table("~/ARG_with_SNP/14880file_ARGanno_absence01_SNPmatrix959090_absence01",tempstr,"_addPhenotype.txt",header = T,sep="\t",quote = '',comment.char = '',check.names = F)
data_arg<-data_all[,c(1,18:ncol(data_all),num)]
delindex<-which(data_arg[,ncol(data_arg)]=="N/A")
if(length(delindex)>0){
  data_arg2<-data_arg[-delindex,]
}else{
  data_arg2<-data_arg
}
delindex2<-which(data_arg2[,ncol(data_arg)]=="I")
if(length(delindex2)>0){
  data_arg3<-data_arg2[-delindex2,]
}else{
  data_arg3<-data_arg2
}
rownames(data_arg3)<-data_arg3[,1]
data_arg3<-data_arg3[,-1]
colnames(data_arg3)[ncol(data_arg3)]<-"Group"

data_arg3_R<-data_arg3[which(data_arg3[,ncol(data_arg3)]=="R"),]
data_arg3_S<-data_arg3[which(data_arg3[,ncol(data_arg3)]=="S"),]
testIndexR<-sample(1:nrow(data_arg3_R),size=(nrow(data_arg3_R)*1/5),replace = F)
testIndexS<-sample(1:nrow(data_arg3_S),size=(nrow(data_arg3_S)*1/5),replace = F)
data_arg3_R_test<-data_arg3_R[testIndexR,]
data_arg3_S_test<-data_arg3_S[testIndexS,]
data_arg3_test<-rbind(data_arg3_R_test,data_arg3_S_test)
write.table(data_arg3_test,paste0("14880file_ARGanno_absence01_SNPmatrix959090_absence01_addPhenotype_",colnames(data_arg)[ncol(data_arg)],"_verify_test.txt"),sep="\t",quote=F,row.names = T,col.names = T)

data_arg3_2<-rbind(data_arg3_R[-testIndexR,],data_arg3_S[-testIndexS,])
data_arg3_2_R<-data_arg3_2[which(data_arg3_2[,ncol(data_arg3_2)]=="R"),]
data_arg3_2_S<-data_arg3_2[which(data_arg3_2[,ncol(data_arg3_2)]=="S"),]
data_arg3_2_train<-rbind(data_arg3_2_R,data_arg3_2_S)
write.table(data_arg3_2_train,paste0("14880file_ARGanno_absence01_SNPmatrix959090_absence01_addPhenotype_",colnames(data_arg)[ncol(data_arg)],"_train.txt"),sep="\t",quote=F,row.names = T,col.names = T)


ptm2<-proc.time() - ptm
print(ptm2)
