library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DESeq2)
library(DOSE)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gplots)
library(pheatmap)

 
#BCF1&CBF1




###Cast
setwd("/public/home/shidetong/projects/lxx/RNA-seq/haplotype/Cast_M/count")
options(stringsAsFactors = FALSE)
control1<-read.table("Speci_Unique_Cast-RNA-1A_M.count",sep = "\t",col.names = c("gene_id","control1"))
control2<-read.table("Speci_Unique_Cast-RNA-3A_M.count",sep = "\t",col.names = c("gene_id","control2"))
treat1<-read.table("Speci_Unique_Cast-RNA-1A_P.count",sep = "\t",col.names = c("gene_id","treat1"))
treat2<-read.table("Speci_Unique_Cast-RNA-3A_P.count",sep = "\t",col.names = c("gene_id","treat2"))


#B6
setwd("/public/home/shidetong/projects/lxx/RNA-seq/haplotype/B6_M/count")
options(stringsAsFactors = FALSE)
control1<-read.table("Speci_Unique_B6-RNA-1A_M.count",sep = "\t",col.names = c("gene_id","control1"))
control2<-read.table("Speci_Unique_B6-RNA-3A_M.count",sep = "\t",col.names = c("gene_id","control2"))
treat1<-read.table("Speci_Unique_B6-RNA-1A_P.count",sep = "\t",col.names = c("gene_id","treat1"))
treat2<-read.table("Speci_Unique_B6-RNA-1A_P.count",sep = "\t",col.names = c("gene_id","treat2"))


####all
setwd('/public/home/shidetong/projects/lxx/RNA-seq/20210324/sort')
options(stringsAsFactors = FALSE)
control1<-read.table("B6-RNA-1A.count",sep = "\t",col.names = c("gene_id","control1"))
control2<-read.table("B6-RNA-3A.count",sep = "\t",col.names = c("gene_id","control2"))
treat1<-read.table("Cast-RNA-1A.count",sep = "\t",col.names = c("gene_id","treat1"))
treat2<-read.table("Cast-RNA-3A.count",sep = "\t",col.names = c("gene_id","treat2"))


##
# raw_count <- merge(merge(control1, control2, by="gene_id"), merge(treat1, treat2, by="gene_id"))
# raw_count_filt <- raw_count[-1:-5,]
# ENSEMBL <- gsub("\\.\\d*", "", raw_count_filt$gene_id) 
# row.names(raw_count_filt) <- ENSEMBL
# condition <- factor(c(rep("control",2),rep("treat",2)), levels = c("control","treat"))
# mycounts  <- raw_count_filt[2:5]
# colData <- data.frame(row.names=colnames(mycounts), condition)
##
raw_count <- merge(merge(control1, control2, by="gene_id"), merge(treat1, treat2, by="gene_id"))
raw_count_filt <- raw_count[-1:-5,]
ENSEMBL <- gsub("\\.\\d*", "", raw_count_filt$gene_id)
#gene_id  <- raw_count_filt %>% .$gene_id
raw_count_filt$G_id  <- ENSEMBL
#row.names(raw_count_filt) <- ENSEMBL
mycounts  <- raw_count_filt[2:6]
mycounts  <- mycounts[,c(5,1:4)]
# colData <- data.frame(row.names=colnames(mycounts), condition)

data  <- mycounts
data$symbol <- mapIds(org.Mm.eg.db,keys = data$G_id,column = 'SYMBOL',keytype = 'ENSEMBL',multiVals='first')
dataunique <- aggregate(.~symbol,data,max)#去重
dataunique[,3:6] <- apply(dataunique[,3:6],2,function(x){as.numeric(x)})
dataunique <- dataunique[rowSums(dataunique[,3:6])>0,]
row.names(dataunique)  <- dataunique$symbol 
mycounts <- dataunique[3:6]


###
colnames(mycounts) <- c('M1','M2','P1','P2')
write.table(mycounts,'/public/home/shidetong/projects/lxx/RNA-seq/haplotype/B6_M/count/count.csv',row.names=  T,sep = '\t',col.names = T,quote=F)

condition <- factor(c(rep("control",2),rep("treat",2)), levels = c("control","treat"))

condition <- factor(c(rep("control",2),rep("treat",2)), levels = c("treat","control"))
colData <- data.frame(row.names=colnames(mycounts), condition)
mycounts  <- mycounts[1:4]
  

dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=colData, 
                              design= ~ condition)
dds = DESeq(dds)
res = results(dds, contrast=c("condition", "control", "treat"))
res = res[order(res$pvalue),]
summary(res)


table(res$padj<0.05)
table(res$padj<0.05,res$log2FoldChange>0.5)
diff_gene_deseq2 <-subset(res, padj < 0.1 & abs(log2FoldChange) > 0.5)
dim(diff_gene_deseq2)


setwd('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif')
write.table(res,'res.csv',row.names = T,sep = '\t')

setwd('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif')
res = read.table('Cast_gene_pvalue.csv',sep = "\t",check.names = F,header = T)

res <- data.frame(symble = row.names(res), res)
cut_off_pvalue = 0.05
cut_off_logFC = 0.5
res$Sig = ifelse(res$padj < cut_off_pvalue & 
                abs(res$log2FoldChange) >= cut_off_logFC, 
                ifelse(res$log2FoldChange> cut_off_logFC ,'Up','Down'),'None')

p = ggplot(res,aes(x =log2FoldChange,y = -log(padj,10),colour = Sig)) +
    geom_point(alpha = 1,size = 3) +
    scale_color_manual(values=c("blue", "grey","red"))+
    geom_vline(xintercept=c(-1,1),lty=6,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_pvalue),
             lty=4,col="black",lwd=0.8) +
    labs(x="log2(Fold Change)",
       y="-log10 (padj)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.8), 
        legend.position="right", 
        legend.title = element_blank()
    )
res$symble <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange)>=0.5,res$symble,"")#

p+geom_text_repel(data = res,aes(x =res$log2FoldChange,y = -log(padj,10),label = symble ),
                  
                  
                  size  = 3,box.padding = unit(0.3, "lines"),
               
                      point.padding = unit(0.3, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE)





#heatmap
data<-read.table("/public/home/shidetong/wedata/test/JRZ/test1.csv",header = T,index = sep=',')
df<-as.matrix(data)


pheatmap(df$R1,
show_rownames = T,
show_colnames = T,
cluster_cols = F,
cluster_rows=T,
filename='test.pdf',#输出文件的名称
fontsize_row=6, #行字体的大小
height=10,  #输出图片的高度
scale = "row",
angle_col=45, #调整字的角度
color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100),
clustering_distance_rows = 'euclidean', 
clustering_method = 'single',
)

 
data<-read.table("/public/home/shidetong/wedata/test/JRZ/test.csv",header = T,sep=',')
df<-as.matrix(data)




p = ggplot(res,aes(x =log2FoldChange,y = -log(padj,10),colour = Sig)) +
    geom_point(alpha = 1,size = 3) +
    geom_text_repel(
    data = subset(res, padj < 0.000000001 & abs(res$log2FoldChange) >= 0.5),
    aes(label = symble),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
    scale_color_manual(values=c("red", "grey","blue"))+
    geom_vline(xintercept=c(-1,1),lty=6,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_pvalue),
             lty=4,col="black",lwd=0.8) +
    labs(x="log2(Fold Change)",
       y="-log10 (padj)",size = 5)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.8), 
        legend.position="right", 
        legend.title = element_blank(),
    axis.title.x =element_text(size=20),
    axis.title.y =element_text(size=20),
    axis.text =element_text(size=20),
    legend.text = element_text(size=15)
    )
setwd('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif/volplot')
ggsave("volCBF1.png", units="in", dpi=300, width=6, height=6, device="png")