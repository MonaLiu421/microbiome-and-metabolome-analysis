library(vegan)
otu <- read.delim('clipboard', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu<-t(otu)
#observed_species <- estimateR(otu)[1, ]
#observed_species<-as.data.frame(observed_species)
#Chao1  <- estimateR(otu)[2, ]
#Chao1<-as.data.frame(Chao1)
alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    
  
  
  Shannon <- sprintf("%0.2f", Shannon)
  Simpson <- sprintf("%0.2f", Simpson)
  result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson)
}
alpha1 <- alpha_diversity (otu)
write.csv(alpha1,"alpha.csv")
alpha<-read.csv("alpha.csv",header = T,row.names = 1)

library(ggpubr)
library(ggplot2)
library(ggsignif)
#obs
my_comparisons_a <- list(c("Normal","Non-return"))
pdf("obs.pdf",width = 4,height = 4)
p<- ggplot(alpha, aes(x=group, y=observed_species,color=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="observed species")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = my_comparisons_a, label = "p.format")+
  theme_bw()+ scale_colour_manual(values=c("#ff0000","#0084ff"))+
  scale_y_continuous(expand = expand_scale(mult = 0.08))+
  theme(legend.position = "top")+theme(rect=element_rect(fill='white'),
                                      plot.margin=unit(rep(0.5,4), 'lines'),
                                      panel.background=element_rect(fill='transparent', color='black'),
                                      panel.border=element_rect(fill='transparent', color='transparent'),
                                      panel.grid=element_blank())
print(p)
dev.off()
#chao1
pdf("chao.pdf",width = 4,height = 4)
p1<- ggplot(alpha, aes(x=group, y=Chao1,color=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="chao")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = my_comparisons_a, label = "p.format")+
  theme_bw()+ scale_colour_manual(values=c("#ff0000","#0084ff"))+
  scale_y_continuous(expand = expand_scale(mult = 0.08))+
  theme(legend.position = "top")+theme(rect=element_rect(fill='white'),
                                      plot.margin=unit(rep(0.5,4), 'lines'),
                                      panel.background=element_rect(fill='transparent', color='black'),
                                      panel.border=element_rect(fill='transparent', color='transparent'),
                                      panel.grid=element_blank())
print(p1)
dev.off()
#ACE
pdf("ACE.pdf",width = 4,height = 4)
p2<-ggplot(alpha, aes(x=group, y=ACE,color=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="Ace")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = my_comparisons_a, label = "p.format")+
  theme_bw()+ scale_colour_manual(values=c("#ff0000","#0084ff"))+
  scale_y_continuous(expand = expand_scale(mult = 0.08))+
  theme(legend.position = "top")+theme(rect=element_rect(fill='white'),
                                       plot.margin=unit(rep(0.5,4), 'lines'),
                                       panel.background=element_rect(fill='transparent', color='black'),
                                       panel.border=element_rect(fill='transparent', color='transparent'),
                                       panel.grid=element_blank())
print(p2)
dev.off()
#shannon
pdf("shannon.pdf",width = 4,height = 4)
p3<-ggplot(alpha, aes(x=group, y=Shannon,color=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="shannon")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = my_comparisons_a, label = "p.format")+
  theme_bw()+ scale_colour_manual(values=c("#ff0000","#0084ff"))+
  scale_y_continuous(expand = expand_scale(mult = 0.08))+
  theme(legend.position = "top")+theme(rect=element_rect(fill='white'),
                                       plot.margin=unit(rep(0.5,4), 'lines'),
                                       panel.background=element_rect(fill='transparent', color='black'),
                                       panel.border=element_rect(fill='transparent', color='transparent'),
                                       panel.grid=element_blank())
print(p3)
dev.off()
#simpson
pdf("simpson.pdf",width = 4,height = 4)
p4<-ggplot(alpha, aes(x=group, y=Simpson,color=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="simpson")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = my_comparisons_a, label = "p.format")+
  theme_bw()+ scale_colour_manual(values=c("#0073c2", "#efc000"))+
  scale_y_continuous(expand = expand_scale(mult = 0.08))+
  theme(legend.position = "top")
print(p4)
dev.off()
library(patchwork)
p5<-(p/p2)|(p1/p3)
ggsave("alpha.diversity.pdf",p5,width = 6,height = 8)
#======beta====
bray_dis <- vegdist(otu, method = 'bray')
bray_dis<-as.matrix(bray_dis)

design = read.csv("group.csv", header=T)

idx = colnames(bray_dis) %in% rownames(design)  

pcoa = cmdscale(bray_dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, design[match(rownames(points), design[,1]), ])#将行名和第一列匹配，并得出匹配的位置，并将其提取出来


adonis_table = adonis(bray_dis~Group, data=design, permutations = 10000) 

adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
adonis_R2 = adonis_table$aov.tab$R2

adonis_pvalue
adonis_R2

library(ggstar)
pdf("bray Fecal PCoA.pdf",width = 6,height = 4)
p6 = ggplot(points, aes(x=x, y=y, color=Group)) +
  stat_ellipse()+
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="bray-curtis Fecal PCoA")+ggtitle(label="bray-curtis Fecal PCoA",subtitle = "p=0.45       R2=0.01")+
       scale_colour_manual(values=c("#ff0000","#0084ff")) +theme_bw()+theme(rect=element_rect(fill='white'),
                                                                           plot.margin=unit(rep(0.5,4), 'lines'),
                                                                           panel.background=element_rect(fill='transparent', color='black'),
                                                                           panel.border=element_rect(fill='transparent', color='transparent'),
                                                                           panel.grid=element_blank())
print(p6)
dev.off()
