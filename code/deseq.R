library(DESeq2)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
 
setwd("/public/home/shidetong/projects/lxx/RNA-seq/20210324/sort")
options(stringsAsFactors = FALSE)
control1<-read.table("B6-RNA-1A.count",sep = "\t",col.names = c("gene_id","control1"))
control2<-read.table("B6-RNA-3A.count",sep = "\t",col.names = c("gene_id","control2"))
treat1<-read.table("Cast-RNA-1A.count",sep = "\t",col.names = c("gene_id","treat1"))
treat2<-read.table("Cast-RNA-3A.count",sep = "\t",col.names = c("gene_id","treat2"))

raw_count <- merge(merge(control1, control2, by="gene_id"), merge(treat1, treat2, by="gene_id"))
raw_count_filt <- raw_count[-1:-5,]
ENSEMBL <- gsub("\\.\\d*", "", raw_count_filt$gene_id)
#gene_id  <- raw_count_filt %>% .$gene_id
raw_count_filt$G_id  <- ENSEMBL
#row.names(raw_count_filt) <- ENSEMBL
mycounts  <- raw_count_filt[2:6]
mycounts  <- mycounts[,c(5,1:4)]


data  <- mycounts
data$symbol <- mapIds(org.Mm.eg.db,keys = data$G_id,column = 'SYMBOL',keytype = 'ENSEMBL',multiVals='first')

dataunique <- aggregate(.~symbol,data,max)#去重

#删除全为空值的行
dataunique[,3:6] <- apply(dataunique[,3:6],2,function(x){as.numeric(x)})
dataunique <- dataunique[rowSums(dataunique[,3:6])>0,]
row.names(dataunique)  <- dataunique$symbol 
mycounts <- dataunique[3:6]

condition <- factor(c(rep("control",2),rep("treat",2)), levels = c("control","treat"))
colData <- data.frame(row.names=colnames(mycounts), condition)
mycounts  <- mycounts[1:4]

dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=colData, 
                              design= ~ condition)
dds = DESeq(dds)
res = results(dds, contrast=c("condition", "control", "treat"))
res = res[order(res$pvalue),]
summary(res)

table(res$padj<0.01)
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 0.5)
dim(diff_gene_deseq2)

sig.gene= diff_gene_deseq2
sig.gene
gene<-rownames(sig.gene)
head(gene)

gene.df<-bitr(gene, fromType = "SYMBOL", 
              toType = c("ENSEMBL","ENTREZID"),
              OrgDb = org.Mm.eg.db)
#GO
ego_all<-enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "ALL",
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
barplot(ego_all,showCategory = 18,title="The GO_ALL enrichment analysis of all DEGs ")+
scale_size(range = c(5,10))+
scale_y_discrete(labels=function(ego_all)str_wrap(ego_all,width=50))
