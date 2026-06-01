data<-read.table("G:/Table_Figure/正文图-所有抗生素-各国沙门耐药率_add_mutation_models.txt",sep="\t",quote = '',header = T,check.names = F)
#countrylist<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Ecuador","Germany","India","Mexico")
countrylist<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Germany","Mexico")
country_data_use1<-data[which(!is.na(match(data$Country,countrylist))),]
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
brewer.pal(12, "Paired") 
ggplot(country_data_use1,aes(x=Antibiotic,y=Resistance_rates,group=Country,color=Country,fill=Country))+
  geom_bar(stat="identity", position=position_dodge(width=.7),width=.5)+
  theme_bw()+
  #theme_minimal()+
  scale_fill_manual(values = c("Australia"="#A6CEE3","Canada"="#1F78B4","China"="#B2DF8A","USA"="#33A02C","UK"="#FB9A99","Brazil"="#FDBF6F","Chile"="#FF7F00","Ecuador"="#CAB2D6","Germany"="#6A3D9A","India"="#FFFF99","Mexico"="#B15928"))+
  scale_color_manual(values = c("Australia"="#A6CEE3","Canada"="#1F78B4","China"="#B2DF8A","USA"="#33A02C","UK"="#FB9A99","Brazil"="#FDBF6F","Chile"="#FF7F00","Ecuador"="#CAB2D6","Germany"="#6A3D9A","India"="#FFFF99","Mexico"="#B15928"))+
  theme(legend.position = "right",plot.title = element_text(hjust=0.5,vjust=0.5), axis.text = element_text(color = "black"), axis.title=element_text(color = "black"))+
  labs(y="Resistant percentage (%)")
ggsave("正文图-各类抗生素9个国家的耐药率-柱状图_add_mutation_models.pdf",width = 10,height = 5,path="G:/Table_Figure")

country_data_use1$Antibiotic<-factor(country_data_use1$Antibiotic,levels=c("AMP","AMS","CFZ","CTX","CAZ","AZM","GEN","CIP","NAL","CT","CHL","TET","SXT"))
ggplot(country_data_use1,aes(x=Country,y=Resistance_rates,group=Antibiotic,color=Antibiotic,fill=Antibiotic))+
  geom_bar(stat="identity", position=position_dodge(width=0.8),width=.2)+
  theme_bw()+
  #theme_minimal()+
  scale_fill_manual(values = c("AMP"="#59A629","AMS"="#FDBF6F","CFZ"="#AE8052","CTX"="#FD8B6F","CAZ"="#A68A29","AZM"="#519FF2","GEN"="#8652F5","CIP"="#51E2F3","NAL"="#20919E","CT"="#F991F4","CHL"="#E56491","TET"="#D1CB4C","SXT"="#87F96E"))+
  scale_color_manual(values = c("AMP"="#59A629","AMS"="#FDBF6F","CFZ"="#AE8052","CTX"="#FD8B6F","CAZ"="#A68A29","AZM"="#519FF2","GEN"="#8652F5","CIP"="#51E2F3","NAL"="#20919E","CT"="#F991F4","CHL"="#E56491","TET"="#D1CB4C","SXT"="#87F96E"))+
  theme(legend.position = "right",plot.title = element_text(hjust=0.5,vjust=0.5), axis.text = element_text(color = "black"), axis.title=element_text(color = "black"))+
  labs(y="Resistant percentage (%)")
ggsave("正文图-各类抗生素9个国家的耐药率-柱状图_add_mutation_models_20251222.pdf",width = 10,height = 5,path="G:/Table_Figure")


##### I级风险等级占比
data<-read.table("G:/Table_Figure/正文图-所有抗生素-各国沙门I风险等级比率_add_mutation_models.txt",sep="\t",quote = '',header = T,check.names = F)
countrylist<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Germany","Mexico")
country_data_use1<-data[which(!is.na(match(data$Country,countrylist))),]
country_data_use2<-country_data_use1
country_data_use2$I_Ratio[which(country_data_use2$Genome_all_number<30)]=0
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
brewer.pal(12, "Paired") 
ggplot(country_data_use2,aes(x=Antibiotic,y=I_Ratio,group=Country,color=Country,fill=Country))+
  #geom_bar(stat="identity", position=position_dodge(width=.7),width=.5)+
  geom_bar(stat="identity", position=position_dodge(width=.8),width=.2)+
  theme_bw()+
  #theme_minimal()+
  scale_fill_manual(values = c("Australia"="#A6CEE3","Canada"="#1F78B4","China"="#B2DF8A","USA"="#33A02C","UK"="#FB9A99","Brazil"="#FDBF6F","Chile"="#FF7F00","Ecuador"="#CAB2D6","Germany"="#6A3D9A","India"="#FFFF99","Mexico"="#B15928"))+
  scale_color_manual(values = c("Australia"="#A6CEE3","Canada"="#1F78B4","China"="#B2DF8A","USA"="#33A02C","UK"="#FB9A99","Brazil"="#FDBF6F","Chile"="#FF7F00","Ecuador"="#CAB2D6","Germany"="#6A3D9A","India"="#FFFF99","Mexico"="#B15928"))+
  theme(legend.position = "right",plot.title = element_text(hjust=0.5,vjust=0.5), axis.text = element_text(color = "black"), axis.title=element_text(color = "black"))+
  labs(y="Risk level I percentage (%)")
ggsave("正文图-各类抗生素9个国家的I风险等级比率-柱状图_add_mutation_models.pdf",width = 10,height = 5,path="G:/Table_Figure")

