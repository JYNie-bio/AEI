setwd("G:/Table_Figure")
library("openxlsx")
data<-read.xlsx("表1-各类抗生素敏感性表型预测模型特征情况_更新CIPNALCT_20260514.xlsx",sheet="图数据_重要特征各预警指标情况 (2)")
#data2<-t(data)
data$Feature2<-paste(data$Antibiotic,data$Feature,sep ="|")
data$Feature_importance<-100*data$Feature_importance

library(ggplot2)
#data_drug<-data[,c("Antibiotic.drug.class","Antibiotic")]
data_drug<-data.frame(Group=c(rep("Drug Class",nrow(data)),rep("Antibiotic",nrow(data))),Text=c(data$Antibiotic.drug.class,data$Antibiotic))
data_importance<-data[,c("Feature2","Feature_importance","featureAnno_percentage","HumanAnno_percentage")]
data_mge<-data[,c("Feature2","MGEAsso_percentage","ISonplasmid_percentage","Cross_phylum_number","Cross_genus_number")]
data_trend<-data[,c("Feature2","Increasing_atleast_oneWHOregion","Increasing_atleast_4WHOregion","Stable_in_allWHOregion","Decreasing_in_atleast_oneWHOregion_and_noincreasing")]
#data_importance_plot<-data.frame(Feature2=rep(data_importance[,1],times=3),Group=rep(colnames(data_importance)[2:4],each=nrow(data_importance)),Value=c(data_importance[,2],data_importance[,3],data_importance[,4]))

##### 特征所属的抗生素、抗生素药物类别柱图 #######
druglist<-unique(data$Antibiotic.drug.class)
antibioticlist<-unique(data$Antibiotic)
data_drug_plot<-matrix(nrow=(length(druglist)+length(antibioticlist)),ncol=3)
colnames(data_drug_plot)<-c("Group","Text","Number")
data_drug_plot[,1]<-c(rep("Drug Class",length(druglist)),rep("Antibiotic",length(antibioticlist)))
data_drug_plot[,2]<-c(druglist,antibioticlist)
for(i in 1:length(druglist)){
  data_drug_plot[i,3]<-length(which(data$Antibiotic.drug.class==data_drug_plot[i,2]))
}
for(i in (1+length(druglist)):nrow(data_drug_plot)){
  data_drug_plot[i,3]<-length(which(data$Antibiotic==data_drug_plot[i,2]))
}
data_drug_plot<-as.data.frame(data_drug_plot)
data_drug_plot$Text<-factor(data_drug_plot$Text,levels = data_drug_plot$Text)
data_drug_plot$Group<-factor(data_drug_plot$Group,levels=c("Drug Class","Antibiotic"))
data_drug_plot$Number<-as.numeric(data_drug_plot$Number)
library("ggfittext")
p0<-ggplot(data_drug_plot,aes(x=Group,y=Number,group=Text,fill=Text))+geom_bar(stat="identity")+theme_classic()+
  geom_bar_text(position="stack",place="center",color="black",label=data_drug_plot$Text)+
  scale_fill_manual(values=c("beta-lactam"="#FCA65A","fluoroquinolone"="#51E2F3","macrolide"="#51A0F3","aminoglycoside"="#8451F3","peptide"="#F890F5","phenicol"="#F890B0","tetracycline"="#F9F26E","sulfonamide-diaminopyrimidine"="#86F96E",
                             "AMP"="#F5AE24","AMS"="#F78B26","CFZ"="#EB8572","CTX"="#F2966C","CAZ"="#D79659","CIP"="#06A5B8","NAL"="#0FB8B0","AZM"="#51A0F3","GEN"="#8451F3","CT"="#F890F5","CHL"="#F890B0","TET"="#F9F26E","SXT"="#86F96E"))
ggsave("图5A-药物类别柱图.pdf",p0,width=6,height=10)

data_drug_plot2<-data_drug_plot[which(data_drug_plot$Group=="Antibiotic"),]
p0<-ggplot(data_drug_plot2,aes(x=Group,y=Number,group=Text,fill=Text))+geom_bar(stat="identity")+theme_classic()+
  geom_bar_text(position="stack",place="center",color="black",label=data_drug_plot2$Text)+
  scale_fill_manual(values=c("AMP"="#FCA65A","AMS"="#FCA65A","CFZ"="#FCA65A","CTX"="#FCA65A","CAZ"="#FCA65A","CIP"="#51E2F3","NAL"="#51E2F3","AZM"="#51A0F3","GEN"="#8451F3","CT"="#F890F5","CHL"="#F890B0","TET"="#F9F26E","SXT"="#86F96E"))
ggsave("图5A-药物类别柱图_only抗生素.pdf",p0,width=2,height=10)


##### 特征重要性、基因组注释率、人类来源注释率等绘图 #######
data_importance$Feature2<-factor(data_importance$Feature2,levels = rev(data$Feature2))
p1<-ggplot(data_importance,aes(x=Feature_importance,y=Feature2))+geom_bar(stat="identity",width = 0.7,col="#638EC6",fill="#638EC6")+theme_test()+
  #scale_y_continuous(limits = c(0, 100))+
  theme(axis.text.y = element_text(angle = 0,vjust=1,hjust=1,color="black",size=10),axis.text.x = element_text(color="black",size=10),axis.title = element_text(color="black",size=10))
#ggsave("图5A-特征重要性柱图.pdf",p1,width=10,height=3.5)
ggsave("图5A-特征重要性柱图.pdf",p1,width=5,height=10)

p2<-ggplot(data_importance,aes(x=featureAnno_percentage,y=Feature2))+geom_bar(stat="identity",width = 0.7,col="#FFB628",fill="#FFB628")+theme_test()+
  #scale_y_continuous(limits = c(0, 100))+
  theme(axis.text.y = element_text(angle = 0,vjust=1,hjust=1,color="black",size=10),axis.text.x = element_text(color="black",size=10),axis.title = element_text(color="black",size=10))
ggsave("图5A-特征基因组注释率柱图.pdf",p2,width=5,height=10)

p3<-ggplot(data_importance,aes(x=HumanAnno_percentage,y=Feature2))+geom_bar(stat="identity",width = 0.7,col="#63C384",fill="#63C384")+theme_test()+
  #scale_y_continuous(limits = c(0, 100))+
  theme(axis.text.y = element_text(angle = 0,vjust=1,hjust=1,color="black",size=10),axis.text.x = element_text(color="black",size=10),axis.title = element_text(color="black",size=10))
ggsave("图5A-特征人类来源注释率柱图.pdf",p3,width=5,height=10)


##### 特征trend of presence 三角形图 #######
library(tidyr)
# 转换为长格式
data_trend_long <- data_trend %>%
  pivot_longer(cols = -Feature2, names_to = "Variable", values_to = "Value")
data_trend_long$Feature2<-factor(data_trend_long$Feature2,levels = rev(data$Feature2))
p4<-ggplot(data_trend_long, aes(x = Variable, y = Feature2, color = Value)) +
  # 添加方框（背景）
  geom_tile(aes(fill = Value), color = "black", width = 0.9, height = 0.9) +
  # 添加点
  geom_point(aes(color = Value), size = 3,shape=17) +
  # 设置颜色
  scale_fill_manual(values = c("Yes" = "white", "No" = "white")) +  # 方框填充色
  scale_color_manual(values = c("Yes" = "#4C3DA2", "No" = "#E2DC69")) +  # 点颜色
  # 调整主题
  theme_minimal() +
  labs(x = "Feature2", y = "Trend of presence")+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45,vjust=1,hjust=1,color="black",size=10),axis.text.y = element_text(color="black",size=10),axis.title = element_text(color="black",size=10))
ggsave("图5A-特征trend情况三角图.pdf",p4,width=7,height=16)

# # 绘制星星图
# library(ggstar)
# ggplot(data_trend_long, aes(x = Variable, y = Feature2, color = Value)) +
#   # 添加方框（背景）
#   geom_tile(aes(fill = Value), color = "black", width = 0.9, height = 0.9) +
#   # # 添加点
#   # geom_point(aes(color = Value), size = 3) +
#   # 添加五角星（默认就是五角星⭐）
#   geom_star(aes(fill = Value,color=Value), size = 3) +  # color=NA 让边框不显示# starshape, 1=五角星，2=六角星，可选其他形状
#   # 设置颜色
#   scale_fill_manual(values = c("Yes" = "white", "No" = "white")) +  # 方框填充色
#   scale_color_manual(values = c("Yes" = "red", "No" = "lightgray")) +  # 点颜色
#   # 调整主题
#   theme_minimal() +
#   labs(x = "Feature2", y = "Trend of presence")+
#   theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45,vjust=1,hjust=1,color="black",size=10),axis.text.y = element_text(color="black",size=10),axis.title = element_text(color="black",size=10))


##### 特征mge情况点图 #######
library(tidyr)
library(dplyr)
library(cowplot)
# 转换为长格式
data_mge_long <- data_mge %>%
  pivot_longer(cols = -Feature2, names_to = "Variable", values_to = "Value")
data_mge_long$Feature2<-factor(data_mge_long$Feature2,levels = rev(data$Feature2))
# 2. 自定义颜色渐变（从浅绿到橙红）
custom_colors <- c("#68A490", "#EAC282", "#D45432")
# 为每列单独绘制点图（隐藏图例，y轴为Feature2）
plot_list <- lapply(unique(data_mge_long$Variable), function(var) {
  df <- data_mge_long %>% filter(Variable == var)
  min_val <- min(df$Value)
  max_val <- max(df$Value)
  p <- ggplot(df, aes(x = 1, y = Feature2)) +
    geom_tile(color = "black", fill = "white", width = 0.9, height = 0.9) +
    geom_point(aes(color = Value), size = 5, shape = 19) +
    scale_color_gradientn(
      colours = custom_colors,
      limits = c(min_val, max_val),
      guide = "none"  # 隐藏子图图例
    ) +
    labs(
      title = paste0(var, "\n[", round(min_val, 1), ", ", round(max_val, 1), "]"),
      x = NULL, y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10, hjust = 1,color="black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5,color="black")
    )
  # 仅第一个图保留y轴标签
  if (var != unique(data_mge_long$Variable)[1]) {
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
})
# 提取所有图例并合并
legends <- lapply(unique(data_mge_long$Variable), function(var) {
  df <- data_mge_long %>% filter(Variable == var)
  min_val <- min(df$Value)
  max_val <- max(df$Value)
  p <- ggplot(df, aes(x = 1, y = 1, color = Value)) +
    geom_point() +
    scale_color_gradientn(
      colours = custom_colors,
      limits = c(min_val, max_val),
      guide = guide_colorbar(title = var, barheight = unit(2, "cm"))
    ) +
    theme_void() +
    theme(legend.position = "right")
  
  get_legend(p)  # 提取图例对象
})

# 拼合主图和图例
library(patchwork)
main_plots <- wrap_plots(plot_list, nrow = 1)  # 纵向排列
all_legends <- wrap_plots(legends, ncol = 1)   # 纵向排列图例
# 最终组合（主图 + 图例）
p5 <- plot_grid(
  main_plots,
  all_legends,
  nrow = 1,
  rel_widths = c(4, 1)  # 主图宽度:图例宽度 = 4:1
)

# 显示图形
print(p5)
ggsave("图5A-特征mge情况点图.pdf",p5,width=8,height=10)



# 2. 绘制热图（分面，每列独立颜色范围）
# ggplot(data_mge_long, aes(x = Variable, y = Feature2, fill = Value)) +
#   # 热图方块
#   geom_tile(color = "black") +
#   # 分面，每列独立
#   facet_wrap(~ Variable, scales = "free_x", nrow = 1) +
#   # 颜色标尺（viridis渐变色，每列独立范围）
#   scale_fill_viridis_c(
#     option = "plasma",
#     guide = guide_colorbar(title = "Value"),
#     na.value = "white"
#   ) +
#   # 调整主题
#   theme_minimal() +
#   theme(
#     axis.text.x = element_blank(),  # 隐藏x轴标签（用分面标题代替）
#     axis.text.y = element_text(size = 10),
#     panel.spacing = unit(0.2, "cm"),
#     strip.text = element_text(face = "bold", size = 10),
#     panel.grid = element_blank()
#   ) +
#   labs(x = NULL, y = "Feature2", fill = "Value")
