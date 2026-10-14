################################################ 有效位点中差异等位数目 ##############################
setwd("~/chewbbaca")
data<-read.table("distance_matrix_all_samples_nie.txt",sep="\t",header=T,row.names = 1,check.names = F)
data <- as.matrix(data)
mode(data) <- "numeric"
data2<-read.table("~/MIC_20230906_14880_addtraintest.txt",sep="\t",header=T,fill=T,quote='')
result<-matrix(nrow=13,ncol=13)
result[,1]<-gsub("_traintest","",colnames(data2)[17:ncol(data2)])
colnames(result)<-c("Antibiotic","traindata_min_max_distance","testdata_min_max_distance","uniquetestdata_min_max_distance","traintestdata_min_max_distance","trainuniquetestdata_min_max_distance","testuniquetestdata_min_max_distance","traindata_quantile_distance","testdata_quantile_distance","uniquetestdata_quantile_distance","traintestdata_quantile_distance","trainuniquetestdata_quantile_distance","testuniquetestdata_quantile_distance")
for(i in 17:ncol(data2)){
  group0<-data2[,c(1,i)]
  colnames(group0)<-c("label","group")
  group2<-as.data.frame(group0[-which(group0$group=="-"),])
  trainlist<-group2$label[which(group2$group=="train")]
  testlist<-group2$label[-which(group2$group=="train")]
  uniquetestlist<-group2$label[which(group2$group=="uniquetest")]
  traindata<-data[trainlist,trainlist]
  testdata<-data[testlist,testlist]
  uniquetestdata<-data[uniquetestlist,uniquetestlist]
  traintestdata<-data[trainlist,testlist]
  trainuniquetestdata<-data[trainlist,uniquetestlist]
  testuniquetestdata<-data[testlist,uniquetestlist]
  # # ==================== 提取距离用于绘制箱线图 ====================
  # dist_test_internal <- testdata[lower.tri(testdata)]
  # dist_uniquetest_internal <- uniquetestdata[lower.tri(uniquetestdata)]
  # dist_train_vs_test <- as.vector(traintestdata)
  # dist_train_vs_uniquetest <- as.vector(trainuniquetestdata)
  # plot_list <- list(
  #   "Test_Internal" = dist_test_internal,
  #   "UniqueTest_Internal" = dist_uniquetest_internal,
  #   "Train_vs_Test" = dist_train_vs_test,
  #   "Train_vs_UniqueTest" = dist_train_vs_uniquetest
  # )
  # df_plot <- stack(plot_list)
  # colnames(df_plot) <- c("Distance", "Group")
  # p <- ggplot(df_plot, aes(x = Group, y = Distance, fill = Group)) +
  #   geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
  #   geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #   scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A",
  #                                "#33A02C")) +
  #   labs(title = paste0("Genomic Distance Distribution - ", colnames(data2)[i]),
  #        y = "cgMLST Distance", x = "") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #         legend.position = "none")
  # ggsave(paste0("Distance_distribution-ownvalid-primary_split-", colnames(data2)[i], ".pdf"),p, width = 9, height = 6)
  # 
  
  trainmin<-min(traindata)
  trainmax<-max(traindata)
  testmin<-min(testdata)
  testmax<-max(testdata)
  uniquetestmin<-min(uniquetestdata)
  uniquetestmax<-max(uniquetestdata)
  traintestmin<-min(traintestdata)
  traintestmax<-max(traintestdata)
  trainuniquetestmin<-min(trainuniquetestdata)
  trainuniquetestmax<-max(trainuniquetestdata)
  testuniquetestmin<-min(testuniquetestdata)
  testuniquetestmax<-max(testuniquetestdata)
  result[i-16,2]<-paste0("[",trainmin,",",trainmax,"]")
  result[i-16,3]<-paste0("[",testmin,",",testmax,"]")
  result[i-16,4]<-paste0("[",uniquetestmin,",",uniquetestmax,"]")
  result[i-16,5]<-paste0("[",traintestmin,",",traintestmax,"]")
  result[i-16,6]<-paste0("[",trainuniquetestmin,",",trainuniquetestmax,"]")
  result[i-16,7]<-paste0("[",testuniquetestmin,",",testuniquetestmax,"]")
  trainquantile<-quantile(traindata, probs = c(0.25, 0.5, 0.75, 0.9))
  testquantile<-quantile(testdata, probs = c(0.25, 0.5, 0.75, 0.9))
  uniquetestquantile<-quantile(uniquetestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  traintestquantile<-quantile(traintestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  trainuniquetestquantile<-quantile(trainuniquetestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  testuniquetestquantile<-quantile(testuniquetestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  result[i-16,8]<-paste0("[", paste(as.numeric(trainquantile), collapse = ", "), "]")
  result[i-16,9]<-paste0("[", paste(as.numeric(testquantile), collapse = ", "), "]")
  result[i-16,10]<-paste0("[", paste(as.numeric(uniquetestquantile), collapse = ", "), "]")
  result[i-16,11]<-paste0("[", paste(as.numeric(traintestquantile), collapse = ", "), "]")
  result[i-16,12]<-paste0("[", paste(as.numeric(trainuniquetestquantile), collapse = ", "), "]")
  result[i-16,13]<-paste0("[", paste(as.numeric(testuniquetestquantile), collapse = ", "), "]")
}
write.table(result,"各抗生素训练集测试集内部及两组间遗传距离范围结果-ownvalid.txt",quote = F,sep="\t",col.names = T,row.names = F)

################################################ 按st和cgmlst重新拆分之后 ##############################
setwd("~/chewbbaca")
library(ggplot2)
data<-read.table("distance_matrix_all_samples_nie.txt",sep="\t",header=T,row.names = 1,check.names = F)
data <- as.matrix(data)
mode(data) <- "numeric"

data2_1<-read.table("重新拆分数据建模/AMP_alldata.txt",sep="\t",header=T,fill=T,quote='')
data2_2<-read.table("重新拆分数据建模/CAZ_alldata.txt",sep="\t",header=T,fill=T,quote='')
data2_3<-read.table("重新拆分数据建模/NAL_alldata.txt",sep="\t",header=T,fill=T,quote='')
data2_12<-cbind(matrix("AMP",nrow=nrow(data2_1)),data2_1[,c(1,ncol(data2_1))])
data2_22<-cbind(matrix("CAZ",nrow=nrow(data2_2)),data2_2[,c(1,ncol(data2_2))])
data2_32<-cbind(matrix("NAL",nrow=nrow(data2_3)),data2_3[,c(1,ncol(data2_3))])
colnames(data2_12)[1]<-"Antibiotic"
colnames(data2_22)[1]<-"Antibiotic"
colnames(data2_32)[1]<-"Antibiotic"
data2<-rbind(rbind(data2_12,data2_22),data2_32)

antibioticlist<-unique(data2$Antibiotic)
result<-matrix(nrow=length(antibioticlist),ncol=13)
result[,1]<-antibioticlist
colnames(result)<-c("Antibiotic","traindata_min_max_distance","testdata_min_max_distance","uniquetestdata_min_max_distance","traintestdata_min_max_distance","trainuniquetestdata_min_max_distance","testuniquetestdata_min_max_distance","traindata_quantile_distance","testdata_quantile_distance","uniquetestdata_quantile_distance","traintestdata_quantile_distance","trainuniquetestdata_quantile_distance","testuniquetestdata_quantile_distance")
plotdata<-matrix()
for(i in 1:nrow(result)){
  group0<-data2[which(data2$Antibiotic==antibioticlist[i]),2:3]
  colnames(group0)<-c("label","group")
  group2<-as.data.frame(group0)
  trainlist<-as.character(group2$label[which(group2$group=="train")])
  testlist<-as.character(group2$label[-which(group2$group=="train")])
  uniquetestlist<-as.character(group2$label[which(group2$group=="uniquetest")])
  traindata<-data[trainlist,trainlist]
  testdata<-data[testlist,testlist]
  uniquetestdata<-data[uniquetestlist,uniquetestlist]
  traintestdata<-data[trainlist,testlist]
  trainuniquetestdata<-data[trainlist,uniquetestlist]
  testuniquetestdata<-data[testlist,uniquetestlist]
  
  # # ==================== 提取距离用于绘制箱线图 ====================
  # dist_test_internal <- testdata[lower.tri(testdata)]
  # dist_uniquetest_internal <- uniquetestdata[lower.tri(uniquetestdata)]
  # dist_train_vs_test <- as.vector(traintestdata)
  # dist_train_vs_uniquetest <- as.vector(trainuniquetestdata)
  # plot_list <- list(
  #   "Test_Internal" = dist_test_internal,
  #   "UniqueTest_Internal" = dist_uniquetest_internal,
  #   "Train_vs_Test" = dist_train_vs_test,
  #   "Train_vs_UniqueTest" = dist_train_vs_uniquetest
  # )
  # df_plot <- stack(plot_list)
  # colnames(df_plot) <- c("Distance", "Group")
  # p <- ggplot(df_plot, aes(x = Group, y = Distance, fill = Group)) +
  #   geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
  #   geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #   scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A",
  #                                "#33A02C")) +
  #   labs(title = paste0("Genomic Distance Distribution - ", antibioticlist[i]),
  #        y = "cgMLST Distance", x = "") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #         legend.position = "none")
  # ggsave(paste0("Distance_distribution-ownvalid-stcgmlst_split-", antibioticlist[i], ".pdf"),p, width = 9, height = 6)
  
  trainmin<-min(traindata)
  trainmax<-max(traindata)
  testmin<-min(testdata)
  testmax<-max(testdata)
  uniquetestmin<-min(uniquetestdata)
  uniquetestmax<-max(uniquetestdata)
  traintestmin<-min(traintestdata)
  traintestmax<-max(traintestdata)
  trainuniquetestmin<-min(trainuniquetestdata)
  trainuniquetestmax<-max(trainuniquetestdata)
  testuniquetestmin<-min(testuniquetestdata)
  testuniquetestmax<-max(testuniquetestdata)
  result[i,2]<-paste0("[",trainmin,",",trainmax,"]")
  result[i,3]<-paste0("[",testmin,",",testmax,"]")
  result[i,4]<-paste0("[",uniquetestmin,",",uniquetestmax,"]")
  result[i,5]<-paste0("[",traintestmin,",",traintestmax,"]")
  result[i,6]<-paste0("[",trainuniquetestmin,",",trainuniquetestmax,"]")
  result[i,7]<-paste0("[",testuniquetestmin,",",testuniquetestmax,"]")
  trainquantile<-quantile(traindata, probs = c(0.25, 0.5, 0.75, 0.9))
  testquantile<-quantile(testdata, probs = c(0.25, 0.5, 0.75, 0.9))
  uniquetestquantile<-quantile(uniquetestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  traintestquantile<-quantile(traintestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  trainuniquetestquantile<-quantile(trainuniquetestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  testuniquetestquantile<-quantile(testuniquetestdata, probs = c(0.25, 0.5, 0.75, 0.9))
  result[i,8]<-paste0("[", paste(as.numeric(trainquantile), collapse = ", "), "]")
  result[i,9]<-paste0("[", paste(as.numeric(testquantile), collapse = ", "), "]")
  result[i,10]<-paste0("[", paste(as.numeric(uniquetestquantile), collapse = ", "), "]")
  result[i,11]<-paste0("[", paste(as.numeric(traintestquantile), collapse = ", "), "]")
  result[i,12]<-paste0("[", paste(as.numeric(trainuniquetestquantile), collapse = ", "), "]")
  result[i,13]<-paste0("[", paste(as.numeric(testuniquetestquantile), collapse = ", "), "]")
}
write.table(result,"各抗生素训练集测试集内部及两组间遗传距离范围结果-ownvalid-重新根据stcgmlst拆分数据后.txt",quote = F,sep="\t",col.names = T,row.names = F)
