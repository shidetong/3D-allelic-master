.libPaths("/public/home/tangy/miniconda3/envs/r413/lib/R/library")
library("DiffBind")
library(ggplot2)
library(parallel)
# setwd('/public/home/shidetong/projects/guolab/now_bow/difbind/JQ1_SW780')
# tamoxifen <- dba(sampleSheet='SW780_JQ1.csv')#读文件，反馈信息列表
# setwd('/public/home/shidetong/projects/guolab/now_bow/difbind/A485_SW780')
# tamoxifen <- dba(sampleSheet='SW780_A485.csv')#读文件，反馈信息列表
# setwd('/public/home/shidetong/projects/guolab/now_bow/difbind/A485_JQ1')
# tamoxifen <- dba(sampleSheet='A485_JQ1.csv')#读文件，反馈信息列表
setwd('/public/home/shidetong/projects/sdt/chip-seq/haplotype/diffbind/number4/M_P')

tamoxifen <- dba(sampleSheet='CBF1BCF1.csv')#读文件，反馈信息列表

#plot(tamoxifen)#一个相关性热图
tamoxifen <- dba.count(tamoxifen)#统计reads数目和反馈FRiP

info <- dba.show(tamoxifen)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
#数据归一化处理
tamoxifen <- dba.normalize(tamoxifen)#
norm <- dba.normalize(tamoxifen, bRetrieve=TRUE)
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,NormLibSize=round(norm$lib.sizes/norm$norm.factors))
#在做差异之前，先构建模型
# tamoxifen <- dba.contrast(tamoxifen,categories=DBA_CONDITION,minMembers = 4)
tamoxifen <- dba.contrast(tamoxifen,categories=DBA_TISSUE,minMembers = 4)


# test <- dba.contrast(tamoxifen,categories=DBA_TISSUE,minMembers = 2)
# test1 <-dba.contrast(tamoxifen,reorderMeta = list(Tissue = 'SW780'),minMembers = 2)

# a <- dba.contrast(tamoxifen,reorderMeta = list(Tissue = 'SW780'),minMembers = 2)

#差异分析
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
# test_a <- dba.analyze(test,method=DBA_ALL_METHODS)
#dba.show(test2_a, bContrasts=TRUE)
# a <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
dba.show(tamoxifen, bContrasts=TRUE)
# dba.show(test_a, bContrasts=TRUE)
#检索差异位点
tamoxifen.DB <- dba.report(tamoxifen)
# tamoxifen.DB <- dba.report(test_a)
sum(tamoxifen.DB$Fold>0)
sum(tamoxifen.DB$Fold<0)

comp1.edgeR <- dba.report(tamoxifen, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(tamoxifen, method=DBA_DESEQ2, contrast = 1, th=1)

out <- as.data.frame(comp1.edgeR)
write.table(out, file="B6_Cast_edgeR.txt", sep="\t", quote=FALSE, col.names = NA)
out <- as.data.frame(comp1.deseq)
write.table(out, file="B6_Cast_deseq2.txt", sep="\t", quote=FALSE, col.names = NA)

out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold","p.value")]
write.table(edge.bed, file="B6_Cast_edgeR_sig.bed", sep="\t", quote=FALSE,row.names=F, col.names=T)


out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold","p.value")]
write.table(deseq.bed, file="B6_Cast_deseq2_sig.bed", sep="\t", quote=FALSE, row.names=F, col.names=T)

#Plot
png(file= 'venn.png')
    
dev.off()
png(file=  'PCA.png')
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_TISSUE)
dev.off()
png(file=  'MA.png')
dba.plotMA(tamoxifen)
dev.off()
png(file= 'volcan.png')
dba.plotVolcano(tamoxifen)
dev.off()
png(file= 'Box.png')
pvals <- dba.plotBox(tamoxifen)
dev.off()
png(file= 'heatmap.png')
corvals <- dba.plotHeatmap(tamoxifen)
dev.off()
png(file=  'dif_heatmap.png')
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE,scale="row", colScheme = hmap)
dev.off()