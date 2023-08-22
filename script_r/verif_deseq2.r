#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

# Get parameters for the test
reads_counts                = args[1]#counts tsv
samples_conditions          = args[2]#samples_conditions tsv

#OPEN TSV FILE
countData=read.csv(reads_counts, header = TRUE, sep='\t')
samples_conditions=read.csv(samples_conditions,header=TRUE, sep = ",")
samples_conditions$SAMPLES=as.character(samples_conditions$SAMPLES)
samples_conditions$CONDITIONS=as.factor(samples_conditions$CONDITIONS)
#samples_conditions=data.frame(samples_conditions)
str(samples_conditions)
## DESeq2 ANALYSIS ON ABUNDANCE TABLE

colData=colnames(countData)
rowData=countData[,1]

rownames(countData)=rowData
countData=countData[,-1]

colData=colnames(countData)

# VECTEUR WITH ORDER CONDITION 

cond=c()

# Parcourir chaque élément de la liste
for (element in colData) {
  # Filtrer le tableau pour les lignes où la colonne1 est égale à l'élément
  similar_lines <- samples_conditions[samples_conditions$SAMPLES == element,]
  #head(similar_lines)
  # Récupérer les valeurs correspondantes de la colonne "valeur"
  cond= c(cond,as.character(similar_lines$CONDITIONS))
  
  
}
colData=data.frame(condition=cond)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)

# PCA
#library(ggfortify)

#mcols(dds)$basepairs <- as.matrix(lenth_genes)
#fpkm <- fpkm(dds)


#pca_des <- rbind(fpkm,type=as.character(DataGroups))
#png("deseq2-pca.png", w=1000, h=1000, pointsize=30)
#autoplot(prcomp(log2((t(fpkm))+1)),data=t(pca_des), colour="type", main="PCA",size=10)+ 
#  theme(plot.title = element_text(face="bold",hjust=0.5,size=50),legend.text=element_text(size=30),
#        legend.title=element_blank(),axis.text = element_text(size=30),
#        axis.title=element_text(size=30))+geom_text(aes(label=colnames(fpkm)),size=10)
#dev.off()

#pdf("deseq2-pca.pdf",10,10)
#autoplot(prcomp(log2((t(fpkm))+1)),data=t(pca_des), colour="type", main="PCA",size=10)+ 
#  theme(plot.title = element_text(face="bold",hjust=0.5,size=50),legend.text=element_text(size=20),
#        legend.title=element_blank(),axis.text = element_text(size=20),
#        axis.title=element_text(size=20))+geom_text(aes(label=colnames(fpkm)),size=7)
#dev.off()

# Run the DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Get differential expression results
res <- results(dds)
res <- res[order(res$padj), ]
## Merge with normalized count data

NormCount<- as.data.frame(counts(dds, normalized=TRUE))
resdata <- merge(as.data.frame(res), as.data.frame(NormCount), by="row.names", sort=FALSE)
#resdata <- merge(as.data.frame(res), as.data.frame(fpkm), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
resdata=resdata[,c(1,3,6:13)]
## Write results
resdata=na.omit(resdata)
resdata[resdata$pvalue==0|resdata$padj==0,c("pvalue","padj")]=0.1e-320
colnames(resdata)[c(2,4)]=c("log2FoldChange","padj")
#write.csv(resdata, file="deseq2-diffexpr-results.tsv", sep="\t")
write.table(resdata, file = "deseq2-diffexpr-results.tsv" , sep = "\t", quote = FALSE, row.names = FALSE)


resdata_filter_pvalue <- resdata[resdata$padj < 0.05, ]
resdata_filter_log <- resdata[abs(resdata$log2FoldChange) > 2, ]
resdata_filter_pvalue_log <- resdata_filter_pvalue[abs(resdata_filter_pvalue$log2FoldChange) > 2, ]

write.table(resdata_filter_pvalue, file = "filter_pvalue.tsv" , sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resdata_filter_log, file = "filter_log.tsv" , sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resdata_filter_pvalue_log, file = "filter_pvalue_log.tsv" , sep = "\t", quote = FALSE, row.names = FALSE)

## Volcano plot
#library(ggrepel)
#resdata$color="F"
#for (i in 1:nrow(resdata)) {
#  if (resdata$padj[i]<0.05) {resdata$color[i]="FDR<0.05"}
#  if (abs(resdata$log2FoldChange[i])>1.5) {resdata$color[i]="|LogFC|>1.5"}
#  if (resdata$padj[i]<0.05 & abs(resdata$log2FoldChange[i])>1.5) {resdata$color[i]="both"}
#}
#resdata=resdata[order(resdata$padj,decreasing = T),]
#png("deseq2-volcanoplot.png", 1200, 1000, pointsize=12)
#print(ggplot(resdata, aes(log2FoldChange, -log10(padj)))+
#        geom_point(aes(color = color)) + 
#        scale_color_manual(values = c("F"="black","FDR<0.05"="red","|LogFC|>1.5"="orange","both"="green"),
#                           breaks=c("F","FDR<0.05", "|LogFC|>1.5","both"), name="",
#                           labels=c("F","FDR<0.05", "|LogFC|>1.5","both"),
#                           limits=c("FDR<0.05", "|LogFC|>1.5","both")) +
#        geom_text_repel(
#          data = subset(resdata, abs(log2FoldChange)>3 & padj <.05),
#          aes(label = Gene), size=2, segment.size=0.2) + ggtitle("Volcano Plot") +
#        theme(plot.title = element_text(hjust = 0.5,face="bold",size=20),
#              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#              panel.background = element_blank(), axis.line = element_line(colour = "black"),
#              axis.text=element_text(size=14),axis.title=element_text(size=14),
#              legend.key = element_rect(colour = NA, fill = NA),
#              legend.title = element_blank(),legend.text=element_text(size=14),
#              legend.background = element_rect(fill="white",
#                                               size=0.25, linetype="solid", 
#                                               colour ="black")))
#dev.off()

#pdf("deseq2-volcanoplot.pdf", 15,10)
#print(ggplot(resdata, aes(log2FoldChange, -log10(padj)))+
#        geom_point(aes(color = color)) + 
#        scale_color_manual(values = c("F"="black","FDR<0.05"="red","|LogFC|>1.5"="orange","both"="green"),
#                           breaks=c("F","FDR<0.05", "|LogFC|>1.5","both"), name="",
#                           labels=c("F","FDR<0.05", "|LogFC|>1.5","both"),
#                           limits=c("FDR<0.05", "|LogFC|>1.5","both")) +
#        geom_text_repel(
#          data = subset(resdata, abs(log2FoldChange)>3 & padj <.05),
#          aes(label = Gene), size=2, segment.size=0.2) + ggtitle("Volcano Plot") +
#        theme(plot.title = element_text(hjust = 0.5,face="bold",size=20),
#              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#              panel.background = element_blank(), axis.line = element_line(colour = "black"),
#              axis.text=element_text(size=14),axis.title=element_text(size=14),
#              legend.key = element_rect(colour = NA, fill = NA),
#              legend.title = element_blank(),legend.text=element_text(size=14),
#              legend.background = element_rect(fill="white",
#                                               size=0.25, linetype="solid", 
#                                               colour ="black")))
#dev.off()


# heatmap topgenes
#library(genefilter)
#library(RColorBrewer)
#library(gplots)
#topVarGenes <- head( order( rowVars( fpkm ), decreasing=TRUE ), 100 )
#png("deseq2-heatmap-topVarGenes.png", w=8, h=9, pointsize=20, res=300, units = "in")
#par(cex.main=0.8)
#heatmap.2( fpkm[ topVarGenes, ], cexCol=0.5, cexRow=0.3, offsetRow=-0.4, offsetCol=-0.4, 
#           scale="row", trace="none", dendrogram="none", main="Top 100 Variance Genes Heatmap",
#           Colv=FALSE, col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), srtCol=30,
#           key.par=list(cex=0.6))
#dev.off()
#
#pdf("deseq2-heatmap-topVarGenes.pdf", w=8, h=8)
#heatmap.2( fpkm[ topVarGenes, ], cexCol=0.9, cexRow=0.5, offsetRow=-0.4, offsetCol=-0.4, 
#           scale="row", trace="none", dendrogram="none", main="Top 100 Variance Genes Heatmap",
#           Colv=FALSE, col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), srtCol=30)
#dev.off()

# Sample distance heatmap
#mycols <- brewer.pal(8, "Dark2")[1:length(unique(DataGroups))]
#sampleDists <- as.matrix(dist(t(assay(rld))))
#png("deseq2-heatmap-samples.png", w=1000, h=1000, pointsize=20)
#heatmap.2(as.matrix(sampleDists), key=F, trace="none",
#          col=colorpanel(100, "black", "white"),
#          ColSideColors=mycols[DataGroups], RowSideColors=mycols[DataGroups],
#          margin=c(10, 10), main="Sample Distance Matrix")
#dev.off()

#pdf("deseq2-heatmap-samples.pdf", w=8, h=8)
#heatmap.2(as.matrix(sampleDists), key=F, trace="none",
#          col=colorpanel(100, "black", "white"),
#          ColSideColors=mycols[DataGroups], RowSideColors=mycols[DataGroups],
#          margin=c(10, 10), main="Sample Distance Matrix")
#dev.off()
