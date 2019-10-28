#!/usr/bin/env Rscript
############################
### Deseq2 analysis pipeline ###
############################


args = commandArgs(trailingOnly=TRUE)

### LOAD PACKAGES ###

library(DESeq2)
library(vsn)
library(ggplot2)
#library(AnnotationDbi)
#library(MSigDB)
#library(gskb)
#library(gage)
#library(ggpubr)
#library(RColorBrewer)
#library(pheatmap)
#library(plyr)


#setwd("~/ferrari/PhD_project/reference_datasets/Ferrari_mESC-iNPC_DMSOvsEPZ_fullEpigenomes_MERGED/downstream_analysis/DE_peaks_analysis_pipe/build_pipe/")

### LOAD COUNT TABLE ###
count=read.csv(args[1], sep='\t', header = F, comment.char = '#')
#count=read.csv("H3K4me3_pipe/H3K4me3_counts.counts", sep='\t', header = F, comment.char = '#')

#out_bed=count[,c("V1","V2","V3")]
#out_bed$ID = paste(out_bed$V1,out_bed$V2,out_bed$V3,sep="_")
#out_bed$score = rep(".",dim(count)[1])
#out_bed$strand = rep(".",dim(count)[1])
#write.table(out_bed, paste(args[3],"peaks_Annotation.bed",sep="/"), sep = "\t", col.names=F, row.names =F, quote=F)

rownames(count) = paste(count$V1,count$V2,count$V3,sep="_")
count$V1 = NULL
count$V2 = NULL
count$V3 = NULL



                            #######################
                            ### mRNA processing ###
                            #######################


### prepare design table ###
des_table = read.csv(args[2],sep="\t")
#des_table = read.csv("info_data.txt",sep="\t")
rownames(des_table) = paste(des_table$mark,des_table$condition,des_table$rep,sep="_")
colnames(count) = rownames(des_table)
all(rownames(des_table) == colnames(count))


### create dds object ###
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = des_table,
                              design = ~ rep + condition)

#hist(log2(rowSums(counts(dds))),n=40)


### PREFILTERING ###
#keep = rowSums(counts(dds)) >= 100
keep = apply(counts(dds), 1, function(x) all(x > 1))
dds = dds[keep,]

file_anno_filt = data.frame(chr = sapply(strsplit(rownames(dds),"_"), `[`, 1),
                            start = sapply(strsplit(rownames(dds),"_"), `[`, 2),
                            end = sapply(strsplit(rownames(dds),"_"), `[`, 3),
                            ID = rownames(dds),
                            score = rep(".",dim(dds)[1]),
                            strande = rep(".",dim(dds)[1]))

write.table(file_anno_filt, paste(args[3],"peaks_Annotation.bed",sep="/"), sep = "\t", col.names=F, row.names =F, quote=F)


### SET REFERENCE LEVELS ###
dds$condition = relevel(dds$condition, ref = as.character(des_table$condition[1]))

### SAMPLE CLUSTERING ### 
#########################

rld <- rlog(dds, blind=FALSE)


### Mean/sd plot ###
pdf(paste(args[3],"Mean_SDplot.pdf",sep="/"), height = 5)
meanSdPlot(assay(rld))
dev.off()


### PCA ###
pdf(paste(args[3],"PCA.pdf",sep="/"), height = 5)
plotPCA(rld, ntop=500, intgroup='condition')
dev.off()


### DIFFERENTIAL EXPRESSION ANALYSIS ###

dds = DESeq(dds)

pdf(paste(args[3],"plotDispEsts.pdf",sep="/"), height = 5)
plotDispEsts(dds)
dev.off()

print(resultsNames(dds))

res = results(dds)
resLFC_normal = lfcShrink(dds, coef=tail(resultsNames(dds),n=1), type="normal")
resLFC_apeglm = lfcShrink(dds, coef=tail(resultsNames(dds),n=1), type="apeglm")


### write results to file for H3K27ac ###
resLFC_normal_ordered = as.data.frame(resLFC_normal[order(resLFC_normal$padj),])
write.table(resLFC_normal_ordered, 
            file = paste(args[3],"resLFC_normal.tsv",sep="/"), 
            quote = F, 
            sep = '\t')

resLFC_apeglm_ordered = as.data.frame(resLFC_apeglm[order(resLFC_apeglm$padj),])
write.table(resLFC_apeglm_ordered, 
            file = paste(args[3],"resLFC_apeglm.tsv",sep="/"), 
            quote = F, 
            sep = '\t')


res_ordered = as.data.frame(res[order(res$padj),])
write.table(res_ordered, 
            file = paste(args[3],"res.tsv",sep="/"), 
            quote = F, 
            sep = '\t')


### VISUALIZATION ###
#####################

### MA-plots ###

pdf(paste(args[3],"MA_res.pdf",sep="/"), height = 5)
DESeq2::plotMA(res, alpha= 0.1, ylim=c(-2,2))
dev.off()

pdf(paste(args[3],"MA_res_normal.pdf",sep="/"), height = 5)
DESeq2::plotMA(resLFC_normal, alpha= 0.1, ylim=c(-2,2))
dev.off()

pdf(paste(args[3],"MA_res_apeglm.pdf",sep="/"), height = 5)
DESeq2::plotMA(resLFC_apeglm, alpha= 0.1, ylim=c(-2,2))
dev.off()


