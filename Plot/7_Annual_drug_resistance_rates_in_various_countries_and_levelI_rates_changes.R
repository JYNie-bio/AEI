setwd("G:/Feature_analysis/predict_use_data_add_mutation_models/")
library(tidyr)
library(ggplot2)

############################################ 十年耐药率折线图 ###########################################
### 计算各个国家的耐药率，创建包含国家名称和耐药率的数据框  
alldata<-read.table("Salmonella_enterica_ownuse14846_publicANI96_173348_Allantibiotic_topfeature_matrix_result_add_year_host_country.txt",sep="\t",quote = '',header = T,comment.char = '',check.names = F,na.strings = T)
alldata<-alldata[-which(alldata$country=="Pacific Ocean"),]
alldata$country[which(alldata$country=="Cote d'Ivoire")]="Ivory Coast"
alldata$country[which(alldata$country=="Korea")]="South Korea"
alldata$country[which(alldata$country=="The United Arab Emirates")]="United Arab Emirates"
alldata$country[which(alldata$country=="U.K.")]="UK"

country_data_all<-matrix(nrow=0,ncol=6)
colnames(country_data_all)=c("Country","Year","Resistance_rates","Genome_resistance_number","Genome_all_number","Antibiotic")

## AZM
#antibiotic="AZM"
for(antibiotic in c("AMP","AMS","AZM","CAZ","CFZ","CHL","CIP","CT","CTX","GEN","NAL","SXT","TET")){
  #index<-which(alldata[,which(colnames(alldata)==antibiotic)]=="I" | alldata[,which(colnames(alldata)==antibiotic)]=="N/A")
  index<-which(alldata[,which(colnames(alldata)==antibiotic)]=="N/A")
  if(length(index)>0){
    data<-alldata[-index,]
  }else{
    data<-alldata
  }
  
  ### 计算各个国家2014-2023年逐年的耐药率
  countrylist<-sort(unique(data$country))
  yearlist<-seq(2014,2023)
  country_data<-matrix(nrow=length(countrylist)*length(yearlist),ncol=5)
  colnames(country_data)<-c("Country","Year","Resistance_rates","Genome_resistance_number","Genome_all_number")
  country_data[,1]=rep(countrylist,each=length(yearlist))
  country_data[,2]=rep(yearlist,times=length(countrylist))
  for(i in 1:nrow(country_data)){
    tempindex=which(data$country==country_data[i,1] & data$year==country_data[i,2])
    if(length(tempindex)>0){
      tempdata=data[tempindex,which(colnames(data)==antibiotic)]
      country_data[i,3]=100*length(which(tempdata=="R"))/length(tempdata)
      country_data[i,4]=length(which(tempdata=="R"))
      country_data[i,5]=length(tempdata)
    }else{
      country_data[i,3]=0
      country_data[i,4]=0
      country_data[i,5]=0
    }
  }
  country_data2<-cbind(country_data,matrix(rep(antibiotic,nrow(country_data)),ncol=1))
  colnames(country_data2)=c("Country","Year","Resistance_rates","Genome_resistance_number","Genome_all_number","Antibiotic")
  #write.table(country_data2,paste0("正文图-",antibiotic,"各国近十年逐年的沙门耐药率_add_mutation_models.txt"),sep="\t",quote = F,col.names = T,row.names = F)
  country_data_all<-rbind(country_data_all,country_data2)
  write.table(country_data_all,"G:/Table_Figure/正文图-所有抗生素-各国近十年逐年的沙门耐药率_add_mutation_models.txt",sep="\t",quote = F,col.names = T,row.names = F)
  
  country_data_use<-as.data.frame(country_data)
  country_data_use$Resistance_rates<-as.numeric(country_data_use$Resistance_rates)
  country_data_use$Genome_all_number<-as.numeric(country_data_use$Genome_all_number)
  country_data_use$Year<-as.character(country_data_use$Year)
  
  countrylist9<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Germany","Mexico")
  country_data_use1<-country_data_use[which(!is.na(match(country_data_use$Country,countrylist9))),]
  country_data_use2<-country_data_use1[which(country_data_use1$Genome_all_number>=30),]
  ggplot(country_data_use2,aes(x=Year,y=Resistance_rates,group=Country,color=Country))+geom_line(size=0.65,na.rm=TRUE)+geom_point(size=2,shape=20,na.rm = TRUE)+
    facet_wrap(~ Country, nrow = 9, ncol = 1, scales="fixed",strip.position = "right") + # m 和 n 需要根据你的实际需求设置
    scale_color_manual(values=c("Australia"="#1B6CA8","Canada"="#3C8D5A","China"="#C77C2C","USA"="#7A4EAB","UK"="#C44E52","Brazil"="#2A9D8F","Chile"="#D95F02","Germany"="#4C78A8","Mexico"="#B279A2"))+
    #scale_x_continuous(breaks = yearlist,limits = c(min(yearlist),max(yearlist)))+
    theme_bw(base_size=10)+
    labs(title=paste0("The percentage of ",antibiotic,"'s resistant changes year by year"))+ylab("Resistant percentage (%)")+theme(panel.grid.major.x=element_line(color="gray88",size=0.28),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(color="black",fill=NA,size=0.5),
                                                                                                                                      strip.background = element_blank(),strip.placement = "outside",strip.text.y.right = element_text(color="black",face="plain",size=9,angle = 0),panel.spacing = unit(0.15, "lines"),
                                                                                                                                      legend.position = "none",plot.title = element_text(hjust=0.5,vjust=0.5,size=12),axis.text = element_text(color = "black",size=8), axis.title=element_text(color = "black",size=10))
  ggsave(paste0("正文图-",antibiotic,"各国沙门耐药率十年变化-9大国_add_mutation_models_only30.pdf"),width = 5,height = 3,path="G:/Table_Figure")
  
  country_data_use2<-country_data_use1
  country_data_use2$Resistance_rates[which(country_data_use1$Genome_all_number==0)]<-NaN
  country_data_use2 <- country_data_use2 %>%
    mutate(
      # 标记：Genome_all_number < 30 记为 grey，其余保留国家名
      color_tag = ifelse(Genome_all_number < 30, "grey_group", as.character(Country))
    )
  
  line_data <- country_data_use2 %>%
    arrange(Country, Year) %>%   # 按国家、年份排序，保证连线顺序正确
    group_by(Country) %>%
    mutate(
      x_end = lead(Year),
      y_end = lead(Resistance_rates),
      seg_color = lead(color_tag) # 线段颜色跟随后一个点的标记
    ) %>%
    ungroup() %>%
    drop_na(x_end, y_end)         # 剔除末尾无下一个点的行
  
  ggplot(country_data_use2, aes(x = Year, y = Resistance_rates)) +
    # 逐段独立线段（AI 可单独选中每一段）
    geom_segment(
      data = line_data,
      aes(x = Year, y = Resistance_rates, xend = x_end, yend = y_end, color = seg_color),
      size = 0.5, na.rm = TRUE
    ) +
    # 散点同步变色
    geom_point(aes(color = color_tag), size = 2, shape = 20, na.rm = TRUE) +
    
    facet_wrap(~ Country, nrow = 9, ncol = 1, strip.position = "right") +
    
    # 配色：原有国家色 + 灰色组
    scale_color_manual(
      values = c(
        "Australia"="#1B6CA8",
        "Canada"="#3C8D5A",
        "China"="#C77C2C",
        "USA"="#7A4EAB",
        "UK"="#C44E52",
        "Brazil"="#2A9D8F",
        "Chile"="#D95F02",
        "Germany"="#4C78A8",
        "Mexico"="#B279A2",
        "grey_group" = "#999999"  # 阈值<30 统一灰色
      )
    ) +
    
    theme_bw(base_size=10) +
    labs(
      title = paste0("The percentage of ", antibiotic, "'s resistant changes year by year"),
      y = "Resistant percentage (%)"
    ) +
    theme(
      panel.grid.major.x = element_line(color="gray88", size=0.28),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="black", fill=NA, size=0.5),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.right = element_text(color="black", face="plain", size=9, angle = 0),
      panel.spacing = unit(0.15, "lines"),
      legend.position = "none",
      plot.title = element_text(hjust=0.5, vjust=0.5, size=12),
      axis.text = element_text(color = "black", size=8),
      axis.title = element_text(color = "black", size=10)
    )
    
  ggsave(paste0("正文图-",antibiotic,"各国沙门耐药率十年变化-9大国_add_mutation_models.pdf"),width = 5,height = 3,path="G:/Table_Figure")
  
}



############################################ 十年第I风险等级基因组数目比例地图 ###########################################
### 计算各个国家的耐药率，创建包含国家名称和耐药率的数据框  
### 计算各个国家的I风险等级的基因组数目占比  
alldata<-read.table("All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel_add_year_country_202605.txt",sep="\t",quote = '',header = T,comment.char = '',check.names = F,na.strings = T)
#alldata<-alldata[-which(alldata$country=="Pacific Ocean"),]
alldata$country[which(alldata$country=="Cote d'Ivoire")]="Ivory Coast"
alldata$country[which(alldata$country=="Korea")]="South Korea"
alldata$country[which(alldata$country=="The United Arab Emirates")]="United Arab Emirates"
alldata$country[which(alldata$country=="U.K.")]="UK"

country_data_all<-matrix(nrow = 0,ncol=6)
colnames(country_data_all)=c("Country","Year","I_ratio","I_genome_number","Genome_all_number","Antibiotic")

## AZM
#antibiotic="AZM"
for(antibiotic in c("AMP","AMS","AZM","CAZ","CFZ","CHL","CIP","CT","CTX","GEN","NAL","SXT","TET")){
  data<-alldata
  ### 计算各个国家2014-2023年逐年的第I等级的比率
  countrylist<-sort(unique(data$country))
  yearlist<-seq(2014,2023)
  country_data<-matrix(nrow=length(countrylist)*length(yearlist),ncol=5)
  colnames(country_data)<-c("Country","Year","I_ratio","I_genome_number","Genome_all_number")
  country_data[,1]=rep(countrylist,each=length(yearlist))
  country_data[,2]=rep(yearlist,times=length(countrylist))
  for(i in 1:nrow(country_data)){
    tempindex=which(data$country==country_data[i,1] & data$year==country_data[i,2])
    if(length(tempindex)>0){
      tempdata=data[tempindex,which(colnames(data)==antibiotic)]
      country_data[i,3]=100*length(which(tempdata=="I"))/length(which(tempdata!="-"))
      country_data[i,4]=length(which(tempdata=="I"))
      country_data[i,5]=length(which(tempdata!="-"))
    }else{
      country_data[i,3]=0
      country_data[i,4]=0
      country_data[i,5]=0
    }
  }
  country_data2<-cbind(country_data,matrix(rep(antibiotic,nrow(country_data)),ncol=1))
  colnames(country_data2)=c("Country","Year","I_ratio","I_genome_number","Genome_all_number","Antibiotic")
  #write.table(country_data2,paste0("正文图-",antibiotic,"各国近十年逐年的沙门I风险等级比率_add_mutation_models.txt"),sep="\t",quote = F,col.names = T,row.names = F)
  country_data_all<-rbind(country_data_all,country_data2)
  write.table(country_data_all,"G:/Table_Figure/正文图-所有抗生素-各国近十年逐年的沙门I风险等级比率_add_mutation_models.txt",sep="\t",quote = F,col.names = T,row.names = F)
  
  country_data_use<-as.data.frame(country_data)
  country_data_use$I_ratio<-as.numeric(country_data_use$I_ratio)
  country_data_use$Genome_all_number<-as.numeric(country_data_use$Genome_all_number)
  country_data_use$Year<-as.character(country_data_use$Year)

  countrylist9<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Germany","Mexico")
  country_data_use1<-country_data_use[which(!is.na(match(country_data_use$Country,countrylist9))),]
  country_data_use2<-country_data_use1
  country_data_use2$I_ratio[which(country_data_use1$Genome_all_number==0)]<-NaN

  country_data_use2 <- country_data_use2 %>%
    mutate(
      # 标记：Genome_all_number < 30 记为 grey，其余保留国家名
      color_tag = ifelse(Genome_all_number < 30, "grey_group", as.character(Country))
    )
  
  line_data <- country_data_use2 %>%
    arrange(Country, Year) %>%   # 按国家、年份排序，保证连线顺序正确
    group_by(Country) %>%
    mutate(
      x_end = lead(Year),
      y_end = lead(I_ratio),
      seg_color = lead(color_tag) # 线段颜色跟随后一个点的标记
    ) %>%
    ungroup() %>%
    drop_na(x_end, y_end)         # 剔除末尾无下一个点的行
  
  ggplot(country_data_use2, aes(x = Year, y = I_ratio)) +
    # 逐段独立线段（AI 可单独选中每一段）
    geom_segment(
      data = line_data,
      aes(x = Year, y = I_ratio, xend = x_end, yend = y_end, color = seg_color),
      size = 0.5, na.rm = TRUE
    ) +
    # 散点同步变色
    geom_point(aes(color = color_tag), size = 2, shape = 20, na.rm = TRUE) +
    
    facet_wrap(~ Country, nrow = 9, ncol = 1, strip.position = "right") +
    
    # 配色：原有国家色 + 灰色组
    scale_color_manual(
      values = c(
        "Australia"="#1B6CA8",
        "Canada"="#3C8D5A",
        "China"="#C77C2C",
        "USA"="#7A4EAB",
        "UK"="#C44E52",
        "Brazil"="#2A9D8F",
        "Chile"="#D95F02",
        "Germany"="#4C78A8",
        "Mexico"="#B279A2",
        "grey_group" = "#999999"  # 阈值<30 统一灰色
      )
    ) +
    
    theme_bw(base_size=10) +
    labs(
      title = paste0("The percentage of ", antibiotic, "'s I risk level changes year by year"),
      y = "Risk level I percentage (%)"
    ) +
    theme(
      panel.grid.major.x = element_line(color="gray88", size=0.28),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="black", fill=NA, size=0.5),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.right = element_text(color="black", face="plain", size=9, angle = 0),
      panel.spacing = unit(0.15, "lines"),
      legend.position = "none",
      plot.title = element_text(hjust=0.5, vjust=0.5, size=12),
      axis.text = element_text(color = "black", size=8),
      axis.title = element_text(color = "black", size=10)
    )
  
  ggsave(paste0("正文图-",antibiotic,"各国沙门I风险等级比率十年变化-9大国_add_mutation_models.pdf"),width = 5,height = 3,path="G:/Table_Figure")
  
  country_data_use2<-country_data_use1[which(country_data_use1$Genome_all_number>=30),]
  ggplot(country_data_use2,aes(x=Year,y=I_ratio,group=Country,color=Country))+geom_line(size=0.65,na.rm=TRUE)+geom_point(size=2,shape=20,na.rm = TRUE)+
    facet_wrap(~ Country, nrow = 9, ncol = 1, scales="fixed",strip.position = "right") + # m 和 n 需要根据你的实际需求设置
    scale_color_manual(values=c("Australia"="#1B6CA8","Canada"="#3C8D5A","China"="#C77C2C","USA"="#7A4EAB","UK"="#C44E52","Brazil"="#2A9D8F","Chile"="#D95F02","Germany"="#4C78A8","Mexico"="#B279A2"))+
    #scale_x_continuous(breaks = yearlist,limits = c(min(yearlist),max(yearlist)))+
    theme_bw(base_size=10)+
    labs(title=paste0("The percentage of ",antibiotic,"'s I risk level changes year by year"))+ylab("Risk level I percentage (%)")+theme(panel.grid.major.x=element_line(color="gray88",size=0.28),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(color="black",fill=NA,size=0.5),
                                                                                                                                         strip.background = element_blank(),strip.placement = "outside",strip.text.y.right = element_text(color="black",face="plain",size=9,angle = 0),panel.spacing = unit(0.15, "lines"),
                                                                                                                                         legend.position = "none",plot.title = element_text(hjust=0.5,vjust=0.5,size=12),axis.text = element_text(color = "black",size=8), axis.title=element_text(color = "black",size=10))
  
  ggsave(paste0("正文图-",antibiotic,"各国沙门I风险等级比率十年变化-9大国_add_mutation_models_only30.pdf"),width = 5,height = 3,path="G:/Table_Figure")
  
}
