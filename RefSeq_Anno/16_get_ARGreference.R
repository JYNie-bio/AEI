setwd("~/RefSeq_Anno")
args=commandArgs(T)
data<-read.table(args[1],header = T,sep="\t",quote = '',comment.char = '',fill=T)
data_1<-data[which(data$Annotation_Source=="ARGcombine"),]
data_2<-data[which(data$Annotation_Source=="RefSeq"),]
data_use1<-data_1[,c(1:8,25:33,42,45,46)]
data_use2<-data_2[,c(1:17,42,45,46)]
colnames(data_use1)<-c("Sample","Gene_id","Gene_contig","Gene_type","Gene_start","Gene_end","Gene_strand","Gene_phase",
                       "alignment_length","mismatches","gap_opens","identity","qcov","scov","evalue","bit_score","subject_id",
                       "Annotation_Source","Use_Gene_Name_Modify","Use_Protein_Name_Modify")
colnames(data_use2)<-c("Sample","Gene_id","Gene_contig","Gene_type","Gene_start","Gene_end","Gene_strand","Gene_phase",
                       "alignment_length","mismatches","gap_opens","identity","qcov","scov","evalue","bit_score","subject_id",
                       "Annotation_Source","Use_Gene_Name_Modify","Use_Protein_Name_Modify")
data_use<-rbind(data_use1,data_use2)

protlist<-unique(data_use$Use_Protein_Name_Modify)
#data_use$bit_score=as.numeric(data_use$bit_score)
#data_use$alignment_length=as.numeric(data_use$alignment_length)
#data_use$qcov=as.numeric(data_use$qcov)
#data_use$scov=as.numeric(data_use$scov)
print(class(data_use$bit_score))
print(class(data_use$alignment_length))
print(class(data_use$qov))
print(class(data_use$sov))
#sort_order<-with(data_use,order(qcov,alignment_length,bit_score,decreasing = T))
sort_order<-with(data_use,order(bit_score,alignment_length,qcov,scov,decreasing = T))
data_use_use<-data_use[sort_order,]
result<-matrix(nrow=length(protlist),ncol=20)
for(i in 1:length(protlist)){
  result[i,]<-as.matrix(data_use_use[which(data_use_use$Use_Protein_Name_Modify==protlist[i])[1],],nrow=1)
}
colnames(result)<-colnames(data_use_use)
write.table(result,args[2],col.names = T,row.names = F,sep="\t",quote=F)
