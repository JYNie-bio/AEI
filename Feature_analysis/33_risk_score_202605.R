args = commandArgs(T)

setwd("~/Feature_analysis/predict_use_data_add_mutation_models")

all_score_all<-matrix(nrow=0,ncol=4)
colnames(all_score_all)<-c("Antibiotic","Feature","Feature_importance","Trend_score")
all_genome_score_all<-matrix(nrow=0,ncol=8)
colnames(all_genome_score_all)<-c("Genome","Antibiotic","Feature","Feature_importance","Trend_score","Feature_host","Feature_lineage_mobility","FeatureMGE_position")

for(antibiotic in args[1]){
  predictdata<-read.table(paste0("Salmonella_enterica_ownuse14846_publicANI96_173348_",antibiotic,"_topfeature_matrix_add_year_host_country_onlyR.txt"),sep="\t",quote='',check.names = F,header = T)
  featurenum=ncol(predictdata)-4
  featurelist=colnames(predictdata)[2:(featurenum+1)]
  
  ############################# 特征频率分布打分和跨地区传播打分 ########################################
  all_score<-matrix(nrow=featurenum,ncol=4)
  colnames(all_score)<-c("Antibiotic","Feature","Feature_importance","Trend_score")
  all_score[,1]<-rep(antibiotic,nrow(all_score))
  all_score[,2]<-featurelist
  
  ### 特征重要性
  importancedata<-read.table(paste0("~/Feature_analysis/use_topfeature/XGBoost_",antibiotic,"_topfeature_finalmodel_feature_importances.txt"),sep="\t",quote='',check.names = F,header = T,row.names = 1)
  all_score[,3]<-importancedata[match(all_score[,2],importancedata$Feature),2]
  
  library(ggplot2)
  ### 特征在基因组中出现频率五年趋势打分
  trendscore<-read.table("Antibiotic_features_CAGRscore_2017-2022_thread5_add_mutation_models.txt",sep="\t",quote='',check.names = F,header = T)
  all_score[,4]<-trendscore[match(all_score[,2],trendscore$Feature),2]
  
  all_score<-as.data.frame(all_score)
  all_score_all<-rbind(all_score_all,all_score)
  
  
  ############################# 公共沙门基因组各个基因组风险打分 ########################################
  predictdata<-read.table(paste0("Salmonella_enterica_ownuse14846_publicANI96_173348_",antibiotic,"_topfeature_matrix_add_mge_onlyR_202605.txt"),sep="\t",quote='',check.names = F,header = T)
  predictdata[which(predictdata$Host=="Companion animal"),'Host']="Animal"
  predictdata[which(predictdata$Host=="livestock"),'Host']="Animal"
  predictdata[which(predictdata$Host=="Livestock"),'Host']="Animal"
  predictdata[which(predictdata$Host=="Poultry"),'Host']="Animal"
  predictdata[which(predictdata$Host=="Wild animal"),'Host']="Animal"
  predictdata[which(predictdata$Host=="Aquatic"),'Host']="Other"
  predictdata[which(predictdata$Host=="Environment"),'Host']="Other"
  predictdata[which(predictdata$Host=="Feed"),'Host']="Other"
  predictdata[which(predictdata$Host=="Water"),'Host']="Other"
  predictdata<-predictdata[-which(predictdata$Host=="Unknown"),]
  sort(unique(predictdata$Host))
  #理论上应该只有Animal、Food、Human、Other
  
  genomefeaturelist<-unique(predictdata[,c(1,5)])
  genomelist<-unique(genomefeaturelist$Genome)
  all_genome_score<-matrix(nrow=nrow(genomefeaturelist),ncol=8)
  colnames(all_genome_score)<-c("Genome","Antibiotic","Feature","Feature_importance","Trend_score","Feature_host","Feature_lineage_mobility","FeatureMGE_position")
  all_genome_score[,1]<-genomefeaturelist$Genome
  all_genome_score[,2]<-rep(antibiotic,nrow(all_genome_score))
  all_genome_score[,3]<-genomefeaturelist$Feature
  for(i in 1:nrow(all_genome_score)){
    if(all_genome_score[i,3]=="-"){
      all_genome_score[i,4]<-0
      all_genome_score[i,5]<-0
      all_genome_score[i,6]<-0
      all_genome_score[i,7]<-0
      all_genome_score[i,8]<-0
    }else{
      ##频率、跨地区打分
      tempindex<-which(all_score_all$Antibiotic==all_genome_score[i,2] & all_score_all$Feature==all_genome_score[i,3])
      all_genome_score[i,4]<-all_score_all$Feature_importance[tempindex]
      all_genome_score[i,5]<-all_score_all$Trend_score[tempindex]
      ##宿主打分
      if(length(grep("Human",predictdata$Host[which(predictdata$Genome==all_genome_score[i,1])[1]]))>0){
        all_genome_score[i,6]=8
      }else if(length(grep("Food",predictdata$Host[which(predictdata$Genome==all_genome_score[i,1])[1]]))>0 | length(grep("Animal",predictdata$Host[which(predictdata$Genome==all_genome_score[i,1])[1]]))>0){
        all_genome_score[i,6]=4
      }else{
        all_genome_score[i,6]=0
      }
      #宿主不是食品动物或人，风险值为0
      #宿主有食品或动物，没有人，风险值为4
      #宿主有人，风险值为8
    
      ##可转移性
      tempindex<-which(predictdata$Genome==all_genome_score[i,1] & predictdata$Feature==all_genome_score[i,3])
      tempdata<-predictdata[tempindex,]
      if(unique(tempdata$IsWith_MGE)=="No"){
        all_genome_score[i,7]=0
        all_genome_score[i,8]=0
      }else if(unique(tempdata$IsWith_MGE)=="Yes"){
        #移动性
        if(length(grep("Phylum",unique(tempdata$CrossedLineage)))>0){
          all_genome_score[i,7]=7
        }else if(length(grep("Class",unique(tempdata$CrossedLineage)))>0){
          all_genome_score[i,7]=6
        }else if(length(grep("Order",unique(tempdata$CrossedLineage)))>0){
          all_genome_score[i,7]=5
        }else if(length(grep("Family",unique(tempdata$CrossedLineage)))>0){
          all_genome_score[i,7]=4
        }else if(length(grep("Genus",unique(tempdata$CrossedLineage)))>0){
          all_genome_score[i,7]=3
	##注意沙门氏菌没有跨种的结果，因为病原菌中没有沙门属下其他种。因此匹配到Species的都是只在沙门内部传播的
        }else if(length(grep("Species",unique(tempdata$CrossedLineage)))>0){
          all_genome_score[i,7]=1
        }else if(unique(tempdata$CrossedLineage)=="-"){
          all_genome_score[i,7]=1
        }else{
          print("Unknown crossed lineage:")
          print(unique(tempdata$CrossedLineage))
        }
        #特征和关联的MGE整体在基因组上位置
        if(length(grep("Plasmid",tempdata$FeatureMGE_genome_position))>0){
          all_genome_score[i,8]=1
        }else{
          all_genome_score[i,8]=0
        }
        
      }else{
        print("ERROR have two or more IsWith_MGE in index:")
        print(tempindex)
      }
      ##没有关联MGEs，风险值为0
      ##关联MGEs但不可跨谱系，风险值为1
      ##最高可跨门纲目科属，风险值分别为7、6、5、4、3
      ##若特征有MGEs，判断特征和关联MGEs在基因组位置，若有在Plasmid，移动风险值加1，否则，移动风险值加0。
    }  
  }
  
  ##合并
  all_genome_score_all<-rbind(all_genome_score_all,all_genome_score)
}

write.table(all_score_all,paste0("All_antibiotic_feature_index_score_",antibiotic,"_ownuse14846_publicANI96_173348_onlyR_202605.txt"),col.names = T,row.names = F,sep="\t",quote = F)
write.table(all_genome_score_all,paste0("All_antibiotic_All_genome_feature_index_score_",antibiotic,"_ownuse14846_publicANI96_173348_onlyR_202605.txt"),col.names = T,row.names = F,sep="\t",quote = F)

################################### 生成各基因组最终各抗生素打分 ########################
all_genome_score_all<-as.data.frame(all_genome_score_all)
all_genome_score_all$Feature_importance<-as.numeric(all_genome_score_all$Feature_importance)
all_genome_score_all$Trend_score<-as.numeric(all_genome_score_all$Trend_score)
all_genome_score_all$Feature_host<-as.numeric(all_genome_score_all$Feature_host)
all_genome_score_all$Feature_lineage_mobility<-as.numeric(all_genome_score_all$Feature_lineage_mobility)
all_genome_score_all$FeatureMGE_position<-as.numeric(all_genome_score_all$FeatureMGE_position)

genomelist<-sort(unique(all_genome_score_all$Genome))
antibioticlist<-sort(unique(all_genome_score_all$Antibiotic))
score<-matrix(nrow=length(genomelist),ncol=(length(antibioticlist)+1))
colnames(score)<-c("Genome",antibioticlist)
score[,1]<-genomelist
for(i in 1:nrow(score)){
  for(j in 2:ncol(score)){
    tempindex<-which(all_genome_score_all$Genome==score[i,1] & all_genome_score_all$Antibiotic==colnames(score)[j])
    tempdata<-all_genome_score_all[tempindex,]
    # 计算每行第5列加第6、7、8列的和
    a <- rowSums(tempdata[,5:8])
    # 计算每行第4列乘以a的和并得到总分数
    score[i,j] <- sum(tempdata[,4] * a)
    if(is.na(score[i,j])){
      print("ERROR in row:")
      print(i)
      print("ERROR in column:")
      print(j)
    }
  }
}
write.table(score,paste0("All_antibiotic_All_genome_risk_score_",antibiotic,"_ownuse14846_publicANI96_173348_onlyR_202605.txt"),col.names = T,row.names = F,sep="\t",quote = F)
  
