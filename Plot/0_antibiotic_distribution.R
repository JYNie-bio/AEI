setwd("G:/data")
data<-read.table("G:/data/MIC_20230906_14880.txt",sep="\t",header = T)
apply(data[,2:ncol(data)],2,table)
lapply(data, table)
sapply(data,table)

result<-matrix(nrow=(ncol(data)-1),ncol=8)
colnames(result)<-c("Antibiotic","NA_num","I_num","R_num","S_num","I_prop","R_prop","S_prop")
result[,1]<-colnames(data)[2:ncol(data)]
for(i in 1:nrow(result)){
  result[i,2]<-length(which(data[,(i+1)]=="N/A"))
  result[i,3]<-length(which(data[,(i+1)]=="I"))
  result[i,4]<-length(which(data[,(i+1)]=="R"))
  result[i,5]<-length(which(data[,(i+1)]=="S"))
  number=as.numeric(result[i,3])+as.numeric(result[i,4])+as.numeric(result[i,5])
  result[i,6]<-100*as.numeric(result[i,3])/number
  result[i,7]<-100*as.numeric(result[i,4])/number
  result[i,8]<-100*as.numeric(result[i,5])/number
}
result2<-result[-which(!is.na(match(result[,1],c("CFX","IPM")))),]
write.table(result2,"G:/Table_Figure/各抗生素使用基因组比例分布.txt",sep="\t",col.names=T,row.names=F,quote=F)
library(openxlsx)
write.xlsx(result2,"G:/Table_Figure/各抗生素使用基因组比例分布.xlsx")