############################################### 附图1B 逐年菌株数目 #############################################
setwd("G:/Table_Figure")
metadata<-read.table("G:/data/Salmonella_MIC_usegenome_20240829_metainfo_addMIC_14880use.txt",sep="\t",quote = '',header = T,comment.char = '')
yearnum<-aggregate(metadata$Use_name,by=list(metadata$year),length)
colnames(yearnum)<-c("Year","Genome_number")
yearnum<-yearnum[-which(yearnum$Year=="未知"),]
yearnum$Year<-as.numeric(yearnum$Year)
library(ggbump)
ggplot(yearnum,aes(x=Year,y=Genome_number))+geom_line(color="#42CEE0",size=1)+
  theme_classic()+geom_area(alpha = 0.6, position = "stack",fill = "#95E8F2", size = 0.4)+
  scale_x_continuous(breaks = yearnum$Year)+labs(y="Genome number")+
  theme(axis.text = element_text(colour = "black",size=12),axis.text.x = element_text(angle = 45,hjust = 1),axis.title = element_text(colour = "black",size=16))
ggsave("附图1B_逐年菌株数目.pdf",width = 6,height = 4)

ggplot(yearnum,aes(x=Year,y=Genome_number))+geom_bump(color="#42CEE0",size=1)+
  theme_classic()+
  scale_x_continuous(breaks = yearnum$Year)+labs(y="Genome number")+
  theme(axis.text = element_text(colour = "black",size=12),axis.text.x = element_text(angle = 45,hjust = 1),axis.title = element_text(colour = "black",size=16))
ggsave("附图1B_逐年菌株数目_2.pdf",width = 6,height = 4)


############################################### 附图1C RS分组数目柱状图 #############################################
setwd("G:/Table_Figure")
metadata<-read.table("G:/data/MIC_20230906_14880.txt",header = T,sep="\t",comment.char = '',quote = '')
MIC_matrix<-matrix(nrow = 45,ncol=3)
MIC_matrix[,1]<-rep(colnames(metadata)[2:16],each=3)
MIC_matrix[,2]<-rep(c("S","I","R"),times=15)
colnames(MIC_matrix)<-c("Antibiotic","MIC_group","Genome_number")
for(j in 2:16){
  MIC_matrix[((j-2)*3+1),3]<-length(which(metadata[,j]=="S"))
  MIC_matrix[((j-2)*3+2),3]<-length(which(metadata[,j]=="I"))
  MIC_matrix[((j-2)*3+3),3]<-length(which(metadata[,j]=="R"))
}
MIC_matrix<-as.data.frame(MIC_matrix)
MIC_matrix$Genome_number<-as.numeric(MIC_matrix$Genome_number)
ggplot(data=MIC_matrix,aes(x=Antibiotic,y=Genome_number,group=MIC_group,fill=MIC_group))+geom_bar(stat="identity",position = position_dodge(),width = 0.8)+theme_classic()+
  labs(y="Genome number")+scale_fill_manual(values = c("R"="#F38979", "S"="#60A7E4","I"="gray"),name = "Group") +
  theme(axis.title = element_text(size=16,color="black"),axis.text = element_text(size=12,color="black"),legend.title = element_text(size=16,color="black"),legend.text = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5,size=15,color="black"))
ggsave("G:/Table_Figure/附图1C_MICgroup分布_IRS.pdf",width=14,height=6)
##去除CFX和IPM
MIC_matrix2<-MIC_matrix[-which(MIC_matrix[,1]=="CFX" | MIC_matrix[,1]=="IPM"),]
ggplot(data=MIC_matrix2,aes(x=Antibiotic,y=Genome_number,group=MIC_group,fill=MIC_group))+geom_bar(stat="identity",position = position_dodge(),width = 0.8)+theme_classic()+
  labs(y="Genome number")+scale_fill_manual(values = c("R"="#F38979", "S"="#60A7E4","I"="gray"),name = "Group") +
  theme(axis.title = element_text(size=16,color="black"),axis.text = element_text(size=12,color="black"),legend.title = element_text(size=16,color="black"),legend.text = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5,size=15,color="black"))
ggsave("G:/Table_Figure/附图1C_MICgroup分布_IRS_use.pdf",width=15,height=6)

