#import data
data<-read.csv("richness-module-input.csv",header = T,check.names = F)
data1<-data[,-c(1,2)]
#data1<-as.data.frame(lapply(data1,as.numeric))
#colnames(data1)<-colnames(data2)
library(magrittr)
library(ggplot2)
library(ggcorrplot)
library(corrplot)
library(dplyr) 
library(caret)
library(Hmisc)
library(reshape2)
library(gplots)
cordata<- round(cor(data1,method = "spearman"),2)
cordata1<-cordata[c(1:3),c(4:59)]#相关系数
melted_cordata <- melt(cordata1)
ref3<- rcorr(as.matrix(data1),type = "spearman")
otu_metabolite_p<-ref3$P
otu_metabolite_p1<-otu_metabolite_p[1:3,4:59]
otu_metabolite_p1_adj = p.adjust (otu_metabolite_p1, method = "BH")
dim (otu_metabolite_p1_adj) = dim (otu_metabolite_p1)
rownames (otu_metabolite_p1_adj) = rownames (otu_metabolite_p1)
colnames (otu_metabolite_p1_adj) = colnames (otu_metabolite_p1)
melted_p_adjust <- melt(otu_metabolite_p1_adj)
melted_cordata$p.adjust<- melted_p_adjust$value
write.csv(melted_cordata,"richness-module-corr.csv")
#visualization
#柱状图
richness<-read.csv("obs-module-fdrless0.05.csv",header = T)
library(ggplot2)
p<-ggplot(data = richness,mapping = aes(x = reorder(module,R2), y = R2, fill = fill)) + 
  geom_bar(stat = 'identity',alpha=0.7)+theme_bw()+theme(axis.text.x=element_text(angle=75,color="black",vjust = 0.95,hjust = 0.95,size=8),
                                                         axis.text.y=element_text(size = 12, colour = "black"),
                                                         axis.title.y=element_text(size = 14))+
  scale_fill_manual(values=c("#0084ff","#ff0000"))+
  theme(legend.position = "top")+theme(rect=element_rect(fill='white'),
                                       plot.margin=unit(rep(0.5,4), 'lines'),
                                       panel.background=element_rect(fill='transparent', color='black'),
                                       panel.border=element_rect(fill='transparent', color='transparent'),
                                       panel.grid=element_blank())+ggtitle("observed-species") +
                                       theme(plot.title = element_text(hjust = 0.5)) 
ggsave("obs-module-asso.pdf",height = 8,width = 8)
