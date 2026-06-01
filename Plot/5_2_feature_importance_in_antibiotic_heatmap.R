# ############################################ 所有特征 ################################
# setwd("G:/Table_Figure")
# data<-read.xlsx("G:/Table_Figure/表1-各类抗生素敏感性表型预测模型特征情况_更新CIPNALCT_20251222.xlsx",sheet = "合并SNP前",startRow = 2)
# df<-data[,c(11,2:10)]
# #df<-data[which(!is.na(match(data$Feature,names(which(table(data$Feature)>1))))),c(11,2:10)]
# druglist<-c("penicillin beta-lactam","cephalosporin (first generation)","cephalosporin (third generation)","macrolide antibiotic","aminoglycoside antibiotic","fluoroquinolone antibiotic","peptide antibiotic","phenicol antibiotic","tetracycline antibiotic","sulfonamide antibiotic, diaminopyrimidine antibiotic")
# 
# df_plot<-df[order(match(df$Antibiotic.drug.class,druglist)),]
# ##挑选存在重要性大于等于0.05的特征绘图
# featurelist<-sort(unique(df_plot$Feature))
# #featurelist<-sort(unique(df_plot$Feature[which(df_plot$Feature.importance>=0.05)]))
# #antibioticlist<-unique(df_plot$Antibiotic)
# antibioticlist<-c("AMP","AMS","CFZ","CTX","CAZ","AZM","GEN","CIP","NAL","CT","CHL","TET","SXT")
# df_plot2<-matrix(nrow=length(featurelist),ncol=length(antibioticlist))
# rownames(df_plot2)<-featurelist
# colnames(df_plot2)<-antibioticlist
# for(i in 1:nrow(df_plot2)){
#   for(j in 1:ncol(df_plot2)){
#     if(length(which(df_plot$Antibiotic==colnames(df_plot2)[j] & df_plot$Feature==rownames(df_plot2)[i]))>0){
#       df_plot2[i,j]=df_plot$Feature.importance[which(df_plot$Antibiotic==colnames(df_plot2)[j] & df_plot$Feature==rownames(df_plot2)[i])]
#     }else{
#       df_plot2[i,j]=0
#     }
#   }
# }
# 
# # 绘制热图
# library(ComplexHeatmap)
# # 定义颜色分模块
# # 定义颜色区间，考虑0值为灰色，然后开始使用蓝色到红色渐变
# library(circlize)
# breaks <- c(0, 0.0006,0.05, 0.1, 0.2,1)
# color_func<-c("#F3F4F5","#E6F5C9","#B3E2CD","#66C2A5","#6DA6E4","#3F7CBE")
# color_gradient <- colorRamp2(breaks = breaks, colors = color_func)# df_plot3, pdf("All_antibiotic_allfeature_importance_20260518_2.pdf",width = 12,height = 40)
# #color_gradient <- colorRamp2(c(min(df_plot2),0.5,max(df_plot2)), c("white", "#E06565","red"))# df_plot2, pdf("All_antibiotic_allfeature_importance_20260518_1.pdf",width = 12,height = 40)
# 
# antibioticlist_order<-c("AMP","AMS","CFZ","CTX","CAZ","AZM","GEN","CIP","NAL","CT","CHL","TET","SXT")
# df_plot3<-df_plot2[,antibioticlist_order]
# for(i in 1:nrow(df_plot3)){
#   for(j in 1:ncol(df_plot3)){
#     if(as.numeric(df_plot3[i,j])<0.0006){
#       df_plot3[i,j]=0
#     }else if(as.numeric(df_plot3[i,j])>=0.0006 & as.numeric(df_plot3[i,j])<0.05){
#       df_plot3[i,j]=0.0006
#     }else if(as.numeric(df_plot3[i,j])>=0.05 & as.numeric(df_plot3[i,j])<0.1){
#       df_plot3[i,j]=0.05
#     }else if(as.numeric(df_plot3[i,j])>=0.1 & as.numeric(df_plot3[i,j])<0.2){
#       df_plot3[i,j]=0.1
#     }else{
#       df_plot3[i,j]=0.2
#     }
#   }
# }
# write.table(df_plot3,"图2A热图数据_202605.txt",sep="\t",quote = F,col.names = T,row.names = T)
# 
# #构建生境的纵向分类注释
# annot_df2_col=data.frame(Antibiotic=antibioticlist_order)
# annot_col2<-c("AMP"="#F3BC51","AMS"="#F3BC51","CFZ"="#F3BC51","CTX"="#F3BC51","CAZ"="#F3BC51","AZM"="#51A0F3","GEN"="#8451F3","CIP"="#51E2F3","NAL"="#51E2F3","CT"="#F890F5","CHL"="#F890B0","TET"="#F9F26E","SXT"="#86F96E")
# ha2 <- rowAnnotation(df = annot_df2_col, col=list(Antibiotic=annot_col2),show_legend = T,annotation_name_rot = 90)
# ha <- HeatmapAnnotation(df = annot_df2_col, col = list(Antibiotic=annot_col2), gp = gpar(col = "black"),show_legend = T)
# 
# pdf("All_antibiotic_allfeature_importance_20260518_2.pdf",width = 12,height = 40)
# Heatmap(df_plot3, name="Feature importance",top_annotation =ha ,#column_title = "Antibiotic",
# #Heatmap(df_plot2, name="Feature importance",top_annotation =ha ,#column_title = "Antibiotic",
#         col = color_gradient,
#         #rect_gp = gpar(col = "black", lwd = 0.1),#背景线颜色
#         column_names_rot = 0,#旋转横坐标标签方向
#         #column_title_side ="top",#更改横轴标题位置，允许的值为“top”或“bottom”
#         row_title_side ="right",#更改纵轴标题位置，允许的值为“left”或“right”
#         row_names_side = "right",column_names_side = "top",
#         #         heatmap_legend_param = list(direction = "horizontal"),
#         column_title_gp = gpar(fontsize = 20, fontface = "bold"), #用于绘制列文本的图形参数
#         row_title_gp = gpar(fontsize = 14, fontface = "bold"),#用于绘制行文本的图形参数
#         #fontface的可能值可以是整数或字符串：1 = plain，2 = bold，3 =斜体，4 =粗体斜体。如果是字符串，则有效值为：“plain”，“bold”，“italic”，“oblique”和“bold.italic”
#         show_row_names = T, show_column_names = T,#是否显示横纵坐标轴标签
#         cluster_rows = T, cluster_columns = F,#是否显示行列的聚类树,关闭聚类s
#         show_row_dend = F,
#         column_dend_height = unit(2, "cm"), row_dend_width = unit(1, "cm"),#更改行列聚类树的高度或宽度
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           # # 在热图单元格内显示数值
#           # if(df_plot2[i,j]>0){
#           #   grid.text(sprintf("%.3f", df_plot2[i, j]), x, y,
#           #             just = "center",          # Make sure text is centered in the cell
#           #             gp = gpar(fontsize = 10))
#           # }
#         }
# )
# dev.off()
# 

############################################ 所有特征-某些特征合并为others ################################
筛选所有多重耐药的特征，其余单药耐药的特征筛选重要性超过5%的，其余的单药耐药特征在各自抗生素内合并成Others。
setwd("G:/Table_Figure")
library("openxlsx")
data<-read.xlsx("G:/Table_Figure/表1-各类抗生素敏感性表型预测模型特征情况_更新CIPNALCT_20251222.xlsx",sheet = "合并SNP前",startRow = 2)
df0<-data[,c(11,2:10)]
df0$Effect<-gsub("pleiotropy","multi-drug resistance",df0$Effect)
df0$Effect<-gsub("Monopower","single-drug resistance",df0$Effect)

df1<-df0[which(df0$Effect=="multi-drug resistance"),]
df2_0<-df0[which(df0$Effect=="single-drug resistance"),]

df2_1<-df2_0[which(as.numeric(df2_0$Feature.importance)>=0.05),]
df2_2<-df2_0[which(as.numeric(df2_0$Feature.importance)<0.05),]
result0<-matrix(ncol=ncol(df2_2),nrow=length(unique(df2_2$Antibiotic)))
result0[,1]<-unique(df2_2$Antibiotic)
result0[,2]<-"-"
result0[,3]<-"-"
result0[,4]<-paste(unique(df2_2$Antibiotic),"Others")
result0[,6]<-"-"
result0[,7]<-"-"
result0[,8]<-"-"
result0[,9]<-"-"
result0[,10]<-"single-drug resistance"

for(i in 1:nrow(result0)){
  temp<-df2_2[which(df2_2$Antibiotic==unique(df2_2$Antibiotic)[i]),]
  result0[i,5]<-sum(as.numeric(temp$Feature.importance))
}
colnames(result0)<-colnames(df2_1)
df2_use<-rbind(df2_1,result0)

df3<-rbind(df1,df2_use)
df3$Feature.importance<-as.numeric(df3$Feature.importance)
#df<-data[which(!is.na(match(data$Feature,names(which(table(data$Feature)>1))))),c(11,2:10)]
druglist<-c("penicillin beta-lactam","cephalosporin (first generation)","cephalosporin (third generation)","macrolide antibiotic","aminoglycoside antibiotic","fluoroquinolone antibiotic","peptide antibiotic","phenicol antibiotic","tetracycline antibiotic","sulfonamide antibiotic, diaminopyrimidine antibiotic")

df_plot<-df3[order(match(df3$Antibiotic.drug.class,druglist)),]
##挑选存在重要性大于等于0.05的特征绘图
featurelist<-sort(unique(df_plot$Feature))
antibioticlist<-c("AMP","AMS","CFZ","CTX","CAZ","AZM","GEN","CIP","NAL","CT","CHL","TET","SXT")
df_plot2<-matrix(nrow=length(featurelist),ncol=length(antibioticlist))
rownames(df_plot2)<-featurelist
colnames(df_plot2)<-antibioticlist
for(i in 1:nrow(df_plot2)){
  for(j in 1:ncol(df_plot2)){
    if(length(which(df_plot$Antibiotic==colnames(df_plot2)[j] & df_plot$Feature==rownames(df_plot2)[i]))>0){
      df_plot2[i,j]=df_plot$Feature.importance[which(df_plot$Antibiotic==colnames(df_plot2)[j] & df_plot$Feature==rownames(df_plot2)[i])]
    }else{
      df_plot2[i,j]=0
    }
  }
}
df_plot2_feature<-unique(df3[,c("Feature","Effect")])
write.table(df_plot2,"图2A热图数据_20260518_mergeothers.txt",sep="\t",quote = F,col.names = T,row.names = T)

otherslist<-c("AMS Others","CFZ Others","CAZ Others","AZM Others","GEN Others","CIP Others","NAL Others","CT Others","CHL Others","TET Others","SXT Others")

df_plot3<-df_plot2[c(sort(setdiff(rownames(df_plot2),otherslist)),otherslist),]

# 绘制热图
library(ComplexHeatmap)
# 定义颜色分模块
# 定义颜色区间，考虑0值为灰色，然后开始使用蓝色到红色渐变
library(circlize)
#color_gradient <- colorRamp2(c(as.numeric(min(df_plot2)),0.4,as.numeric(max(df_plot2))), c("white", "#66C2A5","#3F7CBE"))# df_plot2, pdf("All_antibiotic_allfeature_importance_20260518_1_combineother.pdf",width = 12,height = 40)
max_val <- max(df_plot3, na.rm = TRUE)
breaks <- seq(0, max_val, length.out = 9)
color_gradient <- colorRamp2(
  breaks,
  c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
)

antibioticlist_order<-c("AMP","AMS","CFZ","CTX","CAZ","AZM","GEN","CIP","NAL","CT","CHL","TET","SXT")
# df_plot3<-df_plot2[,antibioticlist_order]
# for(i in 1:nrow(df_plot3)){
#   for(j in 1:ncol(df_plot3)){
#     if(as.numeric(df_plot3[i,j])<0.0006){
#       df_plot3[i,j]=0
#     }else if(as.numeric(df_plot3[i,j])>=0.0006 & as.numeric(df_plot3[i,j])<0.05){
#       df_plot3[i,j]=0.0006
#     }else if(as.numeric(df_plot3[i,j])>=0.05 & as.numeric(df_plot3[i,j])<0.1){
#       df_plot3[i,j]=0.05
#     }else if(as.numeric(df_plot3[i,j])>=0.1 & as.numeric(df_plot3[i,j])<0.2){
#       df_plot3[i,j]=0.1
#     }else{
#       df_plot3[i,j]=0.2
#     }
#   }
# }
# write.table(df_plot3,"图2A热图数据_202605_mergeothers.txt",sep="\t",quote = F,col.names = T,row.names = T)

#构建生境的纵向分类注释
annot_df2_col=data.frame(Antibiotic=antibioticlist_order)
annot_col1<-c("multi-drug resistance"="#A7CDF9","single-drug resistance"="#BFEEAA")
# 构建数据框，确保列名为 'Feature_category'
annot_df_row <- data.frame(Feature_category = df_plot2_feature$Effect,
                           row.names = rownames(df_plot2_feature))
ha2 <- rowAnnotation(df = annot_df_row, col = list(Feature_category = annot_col1), show_legend = TRUE, annotation_name_rot = 0)

#annot_col2<-c("AMP"="#F3BC51","AMS"="#F3BC51","CFZ"="#F3BC51","CTX"="#F3BC51","CAZ"="#F3BC51","AZM"="#51A0F3","GEN"="#8451F3","CIP"="#51E2F3","NAL"="#51E2F3","CT"="#F890F5","CHL"="#F890B0","TET"="#F9F26E","SXT"="#86F96E")
annot_col2<-c("AMP"="#F9DFAD","AMS"="#F9DFAD","CFZ"="#F9DFAD","CTX"="#F9DFAD","CAZ"="#F9DFAD","AZM"="#AED4FD","GEN"="#CEBAFC","CIP"="#AFF1F9","NAL"="#AFF1F9","CT"="#F8D7CA","CHL"="#F8CAE2","TET"="#FFFCC9","SXT"="#C7F8BC")

ha <- HeatmapAnnotation(df = annot_df2_col, col = list(Antibiotic=annot_col2), gp = gpar(col = NA),show_legend = T)

rect_style <- gpar(col = "gray92", lwd = 0.25)

pdf("All_antibiotic_allfeature_importance_20260529_combineother.pdf",width = 12,height = 18)
Heatmap(df_plot3, name="Feature importance",top_annotation =ha ,right_annotation = ha2,#column_title = "Antibiotic",
        #Heatmap(df_plot2, name="Feature importance",top_annotation =ha ,#column_title = "Antibiotic",
        col = color_gradient,
        #rect_gp = gpar(col = "gray88", lwd = 0.1),#背景线颜色
        rect_gp = rect_style,
        column_names_rot = 0,#旋转横坐标标签方向
        #column_title_side ="top",#更改横轴标题位置，允许的值为“top”或“bottom”
        row_title_side ="right",#更改纵轴标题位置，允许的值为“left”或“right”
        row_names_side = "right",column_names_side = "top",
        #         heatmap_legend_param = list(direction = "horizontal"),
        column_title_gp = gpar(fontsize = 20, fontface = "bold"), #用于绘制列文本的图形参数
        row_title_gp = gpar(fontsize = 14, fontface = "bold"),#用于绘制行文本的图形参数
        #fontface的可能值可以是整数或字符串：1 = plain，2 = bold，3 =斜体，4 =粗体斜体。如果是字符串，则有效值为：“plain”，“bold”，“italic”，“oblique”和“bold.italic”
        show_row_names = T, show_column_names = T,#是否显示横纵坐标轴标签
        cluster_rows = F, cluster_columns = F,#是否显示行列的聚类树,关闭聚类s
        show_row_dend = F,
        column_dend_height = unit(2, "cm"), row_dend_width = unit(1, "cm"),#更改行列聚类树的高度或宽度
        cell_fun = function(j, i, x, y, width, height, fill) {
          # 在热图单元格内显示数值
          if(df_plot3[i,j]>0){
            grid.text(sprintf("%.3f", df_plot3[i, j]), x, y,
                      just = "center",          # Make sure text is centered in the cell
                      gp = gpar(fontsize = 10))
          }
        }
)
dev.off()
