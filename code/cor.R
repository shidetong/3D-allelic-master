library(ggplot2)
library(pheatmap)
library(ggpubr)
library(corrplot)
##原始RNAseq
fpkm1<-read.csv("/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_batch2_HiC_2_1/Cast_Matrix/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt",head=F,sep = '\t')
head(fpkm1)[,1:2]
fpkm2<-read.csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_batch2_HiC_4_3/Castmatrix/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt',head = F,sep = '\t')
head(fpkm2)[,1:2]
fpkm<-cbind(fpkm1,fpkm2)
f <- fpkm[,c(2,4)]
names(f)<-c('CBF1_1','CBF1_2')
 
fpkm3<-read.csv("/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_batch2_HiC_2_1/Cast_Matrix/Cooler/Traditional_TADs/Traditional_TADs_DI_40K.txt",head=F,sep = '\t')
head(fpkm3)[,1:2]
fpkm4<-read.csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_batch2_HiC_4_3/Castmatrix/Cooler/Traditional_TADs/Traditional_TADs_DI_40K.txt',head = F,sep = '\t')
head(fpkm4)[,1:2]
fpkm<-cbind(fpkm3,fpkm4)
f <- fpkm[,c(2,4)]
names(f)<-c('CBF1_1','CBF1_2')

# head(fpkm)
# f <- fpkm[,c(2,4)]
# f$Cast1 <- abs(f$V2)
# f$Cast2 <- abs(f$V2.1)
# head(f)
# g <- abs(f[,3])
# h <- abs(f[,4])
# f <- cbind(g,h)
# colnames(f)<-c('BCF1_1','BCF1_2')
# head(f)


ggscatter(f,x = 'CBF1_1', y = 'CBF1_2',color = 'red',size = 3,cor.coef = TRUE,cor.coef.size = 20)+
theme(axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank())+
theme(axis.title.x =element_text(size=60), axis.title.y=element_text(size=60))
setwd('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_cor')
ggsave(filename = 'cor_CBF1_TAD.png',dpi=300)
