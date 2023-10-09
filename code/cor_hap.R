library(ggplot2)
library(pheatmap)
library(ggpubr)
library(corrplot)
##原始RNAseq
fpkm<-read.csv("/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_150bp/depth/alldepth.csv",head=T,sep = '\t')
head(fpkm)[,4:7]
fpkm <- fpkm[rowSums(fpkm[,2:9])>0,]
fpkm <- fpkm[,4:7]


les <-fpkm[fpkm[,1]<10000,]
q <- les[les[,2]<10000,]
w <- q[q[,3]<10000,]
e <- w[w[,4]<10000,]



sample_cor<-cor(fpkm)


ggscatter(fpkm,x = "B6M", y = "B6P")
plot(fpkm)+
theme(axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank())


ggscatter(e,x = "B6P", y = "CastM",color = 'red',size = 3,cor.coef = TRUE,cor.coef.size = 20)+
theme(axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank())+
theme(axis.title.x =element_text(size=60), axis.title.y=element_text(size=60))



pheatmap(sample_cor,cluster_rows = F,
         scale="none",
         cluster_cols = F,
         fontsize_row = 10,
         fontsize_col = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #color = colorRampPalette(c("green3", "white", "blue4"))(100),#换颜色
         angle_col = 45)#修改横轴坐标名倾斜度)


pheatmap(sample_cor,angle_col = 45)

par(oma=c(3,3,3,3)) 
par(mar=c(3,8,8,3))
corrplot(sample_cor,type = "lower", tl.col = "black", tl.srt = 90)
corrplot()

