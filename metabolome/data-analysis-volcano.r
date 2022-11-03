#data normalization and Univariate Analysis were performed by online Metaboanalyst(V5.0) https://www.metaboanalyst.ca/MetaboAnalyst/upload/StatUploadView.xhtml
#volcano plot
library(ggpubr)
library(ggthemes)
data<-read.csv("all-feature-volcano.csv",header = T)
ggscatter(data,x="logFC",y="logP") +theme_base()
# add a column group
data$group ="not-significant"
data$group[which((data$p.ajusted< 0.05)& (data$logFC >1))] ="up-regulated"
data$group[which((data$p.ajusted< 0.05)& (data$logFC < -1))] ="down-regulated"
table(data$group)
#add a column lable
data$lable =""
#order by p value
data<-data[order(data$p.ajusted),]
up <-head(data$X[which(data$group =="up-regulated")],10)
down <-head(data$X[which(data$group =="down-regulated")],10)
#combine up and down
data_top10<-c(as.character(up),as.character(down))
data$lable[match(data_top10,data$X)] <- data_top10
#group
data$group1<-""
data$group1[which(data$group =="not-significant")]<-"not-significant"
data$group1[which(data$group !="not-significant")]<-data$classII
library(ggrepel)
p<-ggscatter(data,x="logFC",y="logP",color = "group",size = 1,
          palette = c("#2f5688","#BBBBBB","#CC0000"),
          font.label = 8,label = data$lable,repel = T,
          xlab = "Log2 Fold-Change(Non-return/Normal)",
          ylab = "-Log10(Adjust p-value)",) +theme_bw() +
          theme(legend.position="none")
p
ggsave("volcano.pdf",p,width = 6,height = 4)


