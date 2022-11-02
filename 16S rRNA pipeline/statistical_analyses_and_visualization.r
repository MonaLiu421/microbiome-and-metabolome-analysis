library(MicrobiotaProcess)
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
#data import
otu<-"stool-filter-table.qza"
rep <- "stool-filter-rep-seqs.qza"
tree <- "stool-slivarooted-tree.qza"
tax <- "taxonomy-deblur-sliva.qza"
sample <- "sample-metadata-w.txt"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeqza=tree)
ps <- ps_dada2
#rarefaction curve 
p_rare <- ggrarecurve(obj=ps_dada2, facetnrow =2 ,factorNames="Description", linesize=0.8,
                      indexNames=c("Observe","Chao1","ACE","Shannon", "Simpson", "J"), 
                      chunks=300) +
  theme(legend.spacing.y=unit(0.02,"cm"),
        legend.text=element_text(size=4))+
  theme_bw() + scale_colour_manual(values=c("#ff0000","#0084ff" ))
ggsave("rarecation.pdf",p_rare,width = 10,height = 8)
#Alpha diversity
library(ggsci)
alphaobj <- get_alphaindex(ps_dada2)
head(as.data.frame(alphaobj))
p_alpha <- ggbox(alphaobj, geom="violin", factorNames="Description",
                 facetnrow=2) + 
  scale_fill_manual(values=alpha(c("#ff0000","#0084ff" ),0.8))+
  theme(strip.background = element_rect(colour=NA, fill="grey"),
        axis.text.x =element_text(angle = 25,hjust = 1,vjust = 1))
p_alpha 
ggsave("alpha.diversity.pdf",p_alpha,width = 10,height = 8)

#pcoa
pcoares <- get_pcoa(obj=ps_dada2, 
                    distmethod="bray", method="hellinger")
pcoaplot <- ggordpoint(obj=pcoares, biplot=FALSE,arrowsize=1.5,
                       speciesannot=FALSE,pc = c(1,2),
                       factorNames=c("Description"),
                       ellipse=T,poinsize =3.5,linesize=0.6,ellipse_alpha=1)+
  scale_fill_manual(values=alpha(c("#ff0000","#0084ff" ),1))+ 
  scale_color_manual(values=c("#ff0000","#0084ff")) 
pcoaplot
ggsave("pcoaplot.pdf",pcoaplot,width = 6,height = 4) 
#adnois test
distme <- get_dist(ps_dada2, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(ps_dada2), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Description <- factor(sampleda$Description)
set.seed(1024)
library(vegan)
adores <- adonis(distme ~ Description, data=sampleda, permutation=10000)
data.frame(adores$aov.tab)

#Bray–Curtis dissimilarity of microbiome
wuyishan<-read.csv("bray-distance-input.csv",header = T)
library(ggpubr)
library(ggplot2)
library(ggsignif)
my_comparisons_a <- list(c("Non-return","Normal"),c("Nor-Non","Normal"),c("Nor-Non","Non-return"))
pdf("bray-distance-wu.pdf",width = 4,height = 4)
wuyishan$group1<-factor(wuyishan$group,levels = c("Non-return","Normal","Nor-Non"))
p<- ggplot(wuyishan, aes(x=group1, y=distance,color=group))+
  geom_boxplot(alpha=0.8, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y="Bray–Curtis dissimilarity")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = my_comparisons_a, label = "p.format")+
  theme_bw()+ scale_color_manual(values=c("#ff0000","#efc000","#0084ff" ))+
  scale_y_continuous(expand = expand_scale(mult = 0.08))+
  theme(legend.position = "top")+ theme(rect=element_rect(fill='white'),
                                        plot.margin=unit(rep(0.5,4), 'lines'),
                                        panel.background=element_rect(fill='transparent', color='black'),
                                        panel.border=element_rect(fill='transparent', color='transparent'),
                                        panel.grid=element_blank())
print(p)
dev.off()
