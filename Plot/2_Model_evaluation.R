setwd("G:/Table_Figure")
library(openxlsx)
library(ggplot2)
data<-read.xlsx("沙门文章用图表-20250328/附表2-各抗生素最终用模型结果.xlsx")
ggplot(data,aes(x=Antibiotic,y=Accuracy))+geom_bar(stat="identity",fill="#8CA7D8")+theme_classic()+
  theme_classic()+
  #theme_hc()+
  theme(axis.text=element_text(size=12,colour = "black"),axis.title = element_text(size=16,colour = "black"),legend.text = element_text(size=12,colour = "black"),legend.title = element_text(size=14,colour = "black"))
ggsave("G:/Table_Figure/图1C_模型各评价指标结果_onlyAccuracy_add_mutation_models.pdf",width = 8,height = 4)
