setwd("G:/Feature_analysis/predict_use_data_add_mutation_models")
data<-read.table("All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_int_status_202605.txt",sep="\t",quote = '',header = T,row.names = 1)
data2<-data
for(j in 1:ncol(data2)){
  data2[,j]<-as.numeric(data2[,j])
}
data3<-data2[,c("AMP","AMS","CAZ","CFZ","CTX","AZM","CT","CIP","NAL","CHL","GEN","SXT","TET")]

## 绘制各类抗生素打分热图--整体一个颜色
library(ComplexHeatmap)
# 定义颜色渐变
library(circlize)
mycol <- colorRamp2(c(min(data3, na.rm = TRUE),max(data3, na.rm = TRUE)), c("white", "red"))
breaks<-c(0,1,5000,10000,13000)
breaks<-c(0,1,4000,7000,10000)
colors<-c("white","#FEEEEE","#FE4848","#FF2929","#FF0202")
mycol <- colorRamp2(breaks, colors)



## 绘制各类抗生素打分热图--分level色块，每个色块一个颜色
library(ComplexHeatmap)
# 定义颜色渐变
library(circlize)
mycol <- colorRamp2(c(min(data3, na.rm = TRUE),max(data3, na.rm = TRUE)), c("white", "red"))
breaks<-c(0,1,5000,10000,13000)
breaks<-c(0,1,4000,7000,10000)
colors<-c("white","#FEEEEE","#FE4848","#FF2929","#FF0202")
mycol <- colorRamp2(breaks, colors)

#构建注释的数据框
annot_df1 = data.frame(level_group=rep("Level I",10))
#为数据框中各行创建对应元素的颜色
annot_row1 = list(level_group=c("Level I"="#E74040","Level II"="#5A99E0","Level III"="#5EC12C"))
#构建水平方向注释
ha1 <- rowAnnotation(df = annot_df1, col = annot_row1, annotation_name_rot = 90,show_legend = T,show_annotation_name = F,
                     annotation_legend_param = list(
                       level_group = list(ncol=1,labels_gp=gpar(fontsize=10),title_gp=gpar(fontsize=10))))
#########
annot_df2 = data.frame(level_group=rep("Level II",8))
annot_row2 = list(level_group=c("Level I"="#E74040","Level II"="#5A99E0","Level III"="#5EC12C"))
ha2 <- rowAnnotation(df = annot_df2, col = annot_row2, annotation_name_rot = 90,show_legend = T,show_annotation_name = F,
                     annotation_legend_param = list(level_group = list(ncol=1,labels_gp=gpar(fontsize=10),title_gp=gpar(fontsize=10))))
                       
#########
annot_df3 = data.frame(level_group=rep("Level III",4))
annot_row3 = list(level_group=c("Level I"="#E74040","Level II"="#5A99E0","Level III"="#5EC12C"))
ha3 <- rowAnnotation(df = annot_df3, col = annot_row3, annotation_name_rot = 90,show_legend = T,show_annotation_name = F,
                     annotation_legend_param = list(level_group = list(ncol=1,labels_gp=gpar(fontsize=10),title_gp=gpar(fontsize=10))))


# 拆分数据
data_level3 <- data3[which(as.numeric(rownames(data3))<4), ]
data_level2 <- data3[which(as.numeric(rownames(data3))>=4 & as.numeric(rownames(data3))<12), ]
data_level1 <- data3[which(as.numeric(rownames(data3))>=12), ]

# 定义各自的颜色
col_level3 <- colorRamp2(
  breaks = c(0, 1, 4000, 7000, 10000),
  colors = c("white", "#F0FFE6", "#99FF66", "#66CC33", "#339900")
)
col_level2 <- colorRamp2(
  breaks = c(0, 1, 4000, 7000, 10000),
  colors = c("white", "#E6F5FF", "#66C2FF", "#3399FF", "#0077CC")
)
col_level1 <- colorRamp2(
  breaks = c(0, 1, 4000, 7000, 10000),
  colors = c("white", "#FEEEEE", "#FE4848", "#FF2929", "#FF0202")
)

# 绘制3个热图
ht_level3 <- Heatmap(
  as.matrix(data_level3), # 排除 Region 列
  col = col_level3,      # 颜色渐变
  left_annotation = ha3,
  heatmap_legend_param = list(title="Genome number"),# 图例及标题
  cluster_rows = F,       # 是否聚类行
  cluster_columns = F,    # 是否聚类列
  show_row_names = T,
  show_column_names = TRUE,
  column_names_rot = 0,#旋转横坐标标签方向
  column_names_side = "top",
  row_names_side ="left",#更改纵轴标题位置，允许的值为“left”或“right”
  row_names_gp = gpar(fontsize = 10),  # 行名字体大小
  column_names_gp = gpar(fontsize = 10), # 列名字体大小
  show_column_dend = F,show_row_dend = F,
  rect_gp = gpar(col = "grey", lwd = 1),#热图每个cell加边框
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 在热图单元格内显示数值
    if(data_level3[i,j]>0){
      grid.text(sprintf("%.f", data_level3[i, j]), x, y,
                just = "center",          # Make sure text is centered in the cell
                gp = gpar(fontsize = 8))
    }
  }
)
ht_level2 <-  Heatmap(
  as.matrix(data_level2), # 排除 Region 列
  col = col_level2,      # 颜色渐变
  left_annotation = ha2,
  heatmap_legend_param = list(title="Genome number"),# 图例及标题
  cluster_rows = F,       # 是否聚类行
  cluster_columns = F,    # 是否聚类列
  show_row_names = T,
  show_column_names = TRUE,
  column_names_rot = 0,#旋转横坐标标签方向
  column_names_side = "top",
  row_names_side ="left",#更改纵轴标题位置，允许的值为“left”或“right”
  row_names_gp = gpar(fontsize = 10),  # 行名字体大小
  column_names_gp = gpar(fontsize = 10), # 列名字体大小
  show_column_dend = F,show_row_dend = F,
  rect_gp = gpar(col = "grey", lwd = 1),#热图每个cell加边框
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 在热图单元格内显示数值
    if(data_level2[i,j]>0){
      grid.text(sprintf("%.f", data_level2[i, j]), x, y,
                just = "center",          # Make sure text is centered in the cell
                gp = gpar(fontsize = 8))
    }
  }
)
ht_level1 <-  Heatmap(
  as.matrix(data_level1), # 排除 Region 列
  col = col_level1,      # 颜色渐变
  left_annotation = ha1,
  heatmap_legend_param = list(title="Genome number"),# 图例及标题
  cluster_rows = F,       # 是否聚类行
  cluster_columns = F,    # 是否聚类列
  show_row_names = T,
  show_column_names = TRUE,
  column_names_rot = 0,#旋转横坐标标签方向
  column_names_side = "top",
  row_names_side ="left",#更改纵轴标题位置，允许的值为“left”或“right”
  row_names_gp = gpar(fontsize = 10),  # 行名字体大小
  column_names_gp = gpar(fontsize = 10), # 列名字体大小
  show_column_dend = F,show_row_dend = F,
  rect_gp = gpar(col = "grey", lwd = 1),#热图每个cell加边框
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 在热图单元格内显示数值
    if(data_level1[i,j]>0){
      grid.text(sprintf("%.f", data_level1[i, j]), x, y,
                just = "center",          # Make sure text is centered in the cell
                gp = gpar(fontsize = 8))
    }
  }
)

# 垂直拼接
combined_ht <- ht_level3 %v% ht_level2 %v% ht_level1
### 输出图形
pdf("G:/Table_Figure/图5A-各类抗生素中各类风险打分菌株数目_热图_add_mutation_models_20250413_202605.pdf",width=7,height = 8)
draw(combined_ht)
dev.off()
# 
# #### 绘制各类抗生素中最高风险菌株数目占比
# data<-read.table("All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel.txt",sep="\t",quote = '',header = T,row.names = 1)
# plotdata<-matrix(nrow=ncol(data),ncol=2)
# colnames(plotdata)<-c("Antibiotic","I_percentage")
# plotdata[,1]<-colnames(data)
# for(i in 1:nrow(plotdata)){
#   plotdata[i,2]<-100*length(which(data[,i]=="I"))/length(which(data[,i]!="-"))
# }
# plotdata<-as.data.frame(plotdata)
# plotdata$I_percentage<-as.numeric(plotdata$I_percentage)
# plotdata$Antibiotic<-factor(plotdata$Antibiotic,levels = rev(c("AMP","AMS","CAZ","CFZ","CTX","AZM","CT","CIP","NAL","CHL","GEN","SXT","TET")))
# library(ggplot2)
# ggplot(plotdata,aes(x=Antibiotic,y=I_percentage))+geom_bar(stat = "identity",color="#82B8ED",fill="#82B8ED")+theme_classic()+
#   labs(y="Percentage of risk level I strains (%)")+
#   geom_text(label = round(plotdata$I_percentage,2), hjust = -0.5, color = "black", size = 4) +  # 在每个柱子上添加数字
#   theme(axis.text = element_text(color = "black"))+coord_flip()
# ggsave("G:/Table_Figure/图5B-各抗生素最高风险菌株数目占比.pdf",width = 9,height = 8)
# ggsave("G:/Table_Figure/图5B-各抗生素最高风险菌株数目占比_add_mutation_models.pdf",width = 9,height = 8)
# 

#### 绘制各类抗生素中level I、II、III菌株数目占比
data<-read.table("All_antibiotic_All_genome_risk_score_ownuse14846_publicANI96_173348_onlyR_risklevel_202605.txt",sep="\t",quote = '',header = T,row.names = 1)
plotdata<-matrix(nrow=ncol(data)*3,ncol=3)
colnames(plotdata)<-c("Antibiotic","Level_group","Percentage")
plotdata[,1]<-rep(colnames(data),each=3)
plotdata[,2]<-rep(c("Level I","Level II","Level III"),times=ncol(data))
for(i in 1:nrow(plotdata)){
  plotdata[i,3]<-100*length(which(data[,which(colnames(data)==plotdata[i,1])]==gsub("Level ","",plotdata[i,2])))/length(which(data[,which(colnames(data)==plotdata[i,1])]!="-"))
}
plotdata<-as.data.frame(plotdata)
plotdata$Percentage<-as.numeric(plotdata$Percentage)
plotdata$Antibiotic<-factor(plotdata$Antibiotic,levels = rev(c("AMP","AMS","CAZ","CFZ","CTX","AZM","CT","CIP","NAL","CHL","GEN","SXT","TET")))
library(ggplot2)
ggplot(plotdata,aes(x=Antibiotic,y=Percentage,fill=Level_group,col=Level_group))+geom_bar(stat = "identity")+theme_classic()+
  labs(y="Percentage of risk level strains (%)")+
  scale_fill_manual(values = c("Level I"="#E74040","Level II"="#5A99E0","Level III"="#5EC12C"))+
  scale_color_manual(values = c("Level I"="#E74040","Level II"="#5A99E0","Level III"="#5EC12C"))+
  #geom_text(label = round(plotdata$Percentage,2), hjust = -0.5, color = "black", size = 4) +  # 在每个柱子上添加数字
  theme(axis.text = element_text(color = "black"))+coord_flip()
ggsave("G:/Table_Figure/图5B-各抗生素各类风险菌株数目占比_add_mutation_models_20260514.pdf",width = 6,height = 8)
