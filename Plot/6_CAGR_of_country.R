setwd("G:/Feature_analysis/predict_use_data_add_mutation_models/")

############################################ 十年耐药率折线图 ###########################################
### 计算各个国家的耐药率，创建包含国家名称和耐药率的数据框  
alldata<-read.table("Salmonella_enterica_ownuse14846_publicANI96_173348_Allantibiotic_topfeature_matrix_result_add_year_host_country.txt",sep="\t",quote = '',header = T,comment.char = '',check.names = F,na.strings = T)
alldata<-alldata[-which(alldata$country=="Pacific Ocean"),]
alldata$country[which(alldata$country=="Cote d'Ivoire")]="Ivory Coast"
alldata$country[which(alldata$country=="Korea")]="South Korea"
alldata$country[which(alldata$country=="The United Arab Emirates")]="United Arab Emirates"
alldata$country[which(alldata$country=="U.K.")]="UK"

CAGR_data_all<-matrix(nrow=0,ncol=3)
colnames(CAGR_data_all)=c("Antibiotic","Country","Resistant_ratio_CAGR")

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
  
  ### 计算各个国家2014-2022年逐年的耐药率
  countrylist<-sort(unique(data$country))
  yearlist<-seq(2014,2022)
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
  
  ##### 计算想要的11个国家的耐药率的CAGR值
  country_data<-as.data.frame(country_data)
  #countrylist5<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Ecuador","Germany","India","Mexico")
  countrylist5<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Germany","Mexico")
  country_data_use1<-country_data[which(!is.na(match(country_data$Country,countrylist5))),]
  country_data_use2<-country_data_use1[which(as.numeric(country_data_use1$Genome_all_number)>=30),]
  countrylist=sort(unique(country_data_use2$Country))
  CAGR_data<-matrix(nrow=length(countrylist),ncol=3)
  colnames(CAGR_data)<-c("Antibiotic","Country","Resistant_ratio_CAGR")
  CAGR_data[,1]<-rep(antibiotic,nrow(CAGR_data))
  CAGR_data[,2]<-countrylist
  for(i in 1:nrow(CAGR_data)){
    tempdata2<-country_data_use2[which(country_data_use2$Country==CAGR_data[i,2]),]
    #if(tempdata2[1,3]==0){tempdata2[1,3]==0.0000000001}
    CAGR_data[i,3]<-100*(((as.numeric(tempdata2[nrow(tempdata2),3])/as.numeric(tempdata2[1,3]))^(1/(as.numeric(tempdata2[nrow(tempdata2),2])-as.numeric(tempdata2[1,2]))))-1)
  }
  CAGR_data_all<-rbind(CAGR_data_all,CAGR_data)
}
#write.table(CAGR_data_all,paste0("G:/Table_Figure/正文图-所有抗生素-11大国近十年逐年耐药率的CAGR.txt"),sep="\t",quote = F,col.names = T,row.names = F)
write.table(CAGR_data_all,paste0("G:/Table_Figure/正文图-所有抗生素-9大国近十年逐年耐药率的CAGR_add_mutation_models.txt"),sep="\t",quote = F,col.names = T,row.names = F)


############################################ 十年第I风险等级基因组数目比例地图 ###########################################
### 计算各个国家的耐药率，创建包含国家名称和耐药率的数据框  
### 计算各个国家的I风险等级的基因组数目占比  
alldata<-read.table("All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel_add_year_country_202605.txt",sep="\t",quote = '',header = T,comment.char = '',check.names = F,na.strings = T)
#alldata<-alldata[-which(alldata$country=="Pacific Ocean"),]
alldata$country[which(alldata$country=="Cote d'Ivoire")]="Ivory Coast"
alldata$country[which(alldata$country=="Korea")]="South Korea"
alldata$country[which(alldata$country=="The United Arab Emirates")]="United Arab Emirates"
alldata$country[which(alldata$country=="U.K.")]="UK"

CAGR_data_all<-matrix(nrow = 0,ncol=3)
colnames(CAGR_data_all)=c("Antibiotic","Country","I_ratio_CAGR")

## AZM
#antibiotic="AZM"
for(antibiotic in c("AMP","AMS","AZM","CAZ","CFZ","CHL","CIP","CT","CTX","GEN","NAL","SXT","TET")){
  data<-alldata
  ### 计算各个国家2014-2022年逐年的第I等级的比率
  countrylist<-sort(unique(data$country))
  yearlist<-seq(2014,2022)
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
  
  ##### 计算想要的9个国家的I风险占比的CAGR值
  country_data<-as.data.frame(country_data)
  countrylist5<-c("Australia","Canada","China","USA","UK","Brazil","Chile","Germany","Mexico")
  country_data_use1<-country_data[which(!is.na(match(country_data$Country,countrylist5))),]
  country_data_use2<-country_data_use1[which(as.numeric(country_data_use1$Genome_all_number)>=30),]
  countrylist=sort(unique(country_data_use2$Country))
  CAGR_data<-matrix(nrow=length(countrylist),ncol=3)
  colnames(CAGR_data)<-c("Antibiotic","Country","I_ratio_CAGR")
  CAGR_data[,1]<-rep(antibiotic,nrow(CAGR_data))
  CAGR_data[,2]<-countrylist
  for(i in 1:nrow(CAGR_data)){
    tempdata2<-country_data_use2[which(country_data_use2$Country==CAGR_data[i,2]),]
    #if(tempdata2[1,3]==0){tempdata2[1,3]==0.0000000001}
    CAGR_data[i,3]<-100*(((as.numeric(tempdata2[nrow(tempdata2),3])/as.numeric(tempdata2[1,3]))^(1/(as.numeric(tempdata2[nrow(tempdata2),2])-as.numeric(tempdata2[1,2]))))-1)
  }
  CAGR_data_all<-rbind(CAGR_data_all,CAGR_data)

}
write.table(CAGR_data_all,paste0("G:/Table_Figure/正文图-所有抗生素-9大国I风险等级比率十年变化的CAGR_add_mutation_models.txt"),sep="\t",quote = F,col.names = T,row.names = F)
