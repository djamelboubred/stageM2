#######################################################################
# THIS SCRIPT TAKE IN INPUT A MATRIX OF ABUNDANCES.
#
# USE DESeq2 libraryfor Differential expression Analysis.
#
#######################################################################

library("data.table")
library("foreach")
library("doParallel")
library("DESeq2")

args <- commandArgs(TRUE)

#gene_id,RC1,RC10,RC2,RC3,RC5,RC6,RC7,RC8,RC9,RS1,RS10,RS2,RS3,RS4,RS5,RS6,RS7,RS8
#Oglab_000003|Oglab_000003,0,9,4,0,8,8,0,0,0,5,4,0,6,6,5,0,10,15
#Oglab_000004|Oglab_000004,0,0,1,0,0,2,0,3,0,3,3,0,0,0,2,2,1,2



# Get parameters for the test
reads_counts                = args[1]#counts tsv
samples_conditions          = args[2]#samples_conditions tsv
pvalue_threshold            = args[3]#pvalue_threshold
log2fc_threshold            = args[4]#log2fc_threshold
conditionA                  = args[5]#condition A name
conditionB                  = args[6]#condition B name
output                      =args[7]#output name

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

str(cond)
colData=data.frame(condition=cond)

#RUN DESeq2

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData,
                              design = ~ condition)

dds <- estimateSizeFactors(dds)

dds <- estimateDispersionsGeneEst(dds)

dds <-  tryCatch(
  estimateDispersionsMAP(estimateDispersionsFit(dds)),
  error=function(e){
    cat("Error during estimateDispersionsFit, probably a '2-order-of-magnitude' message, trying suggested alterative method\n")
    dispersions(dds) <- mcols(dds)$dispGeneEst
    dds
  })

dds <- nbinomWaldTest(dds)

#print(resultsNames(dds))
resDESeq2 <- results(dds, pAdjustMethod = "BH", contrast = c("condition",conditionB,conditionA))

NormCount<- as.data.frame(counts(dds, normalized=TRUE))
#dds_res <- results(dds)



#head(NormCount)
#table(dds_NA_res$padj<=pvalue_threshold)

#print(dds)
#str(NormCount)



#str(resDESeq2)

#head(resDESeq2)

# WRITE A TSV WITH THIS FORMAT FOR THE CURRENT CHUNK
# Kmer_ID
# meanA
# meanB
# log2FC
# NormCount
diff_counts=data.frame(ID=rownames(resDESeq2),
                pvalue=resDESeq2$padj,
                meanA=rowMeans(NormCount[,which(samples_conditions$CONDITIONS == conditionA)]),
                meanB=rowMeans(NormCount[,which(samples_conditions$CONDITIONS == conditionB)]),
                log2FC=resDESeq2$log2FoldChange,
                NormCount)

cat("Le nombre de tags total :",nrow(diff_counts),"\n")

diff_counts=diff_counts[which(diff_counts$pvalue<=pvalue_threshold),]# && which(abs(diff_counts$log2FC)>=log2fc_threshold),]
cat("Le nombre de tags selectionné avec le seuil pour la pvalue:",nrow(diff_counts),"\n")
diff_counts=diff_counts[which(abs(diff_counts$log2FC)>=log2fc_threshold),]
cat("Le nombre de tags selectionné avec le seuil pour log2FC:",nrow(diff_counts),"\n")



#diff_counts=diff_counts[which(diff_counts$pvalue)<=pvalue_threshold && which(abs(diff_counts$log2FC)>=log2fc_threshold),]
#diff_counts=diff_counts[diff_counts[which(diff_counts$pvalue)<=pvalue_threshold && which(abs(diff_counts$log2FC)>=log2fc_threshold)],]
#diff_counts=data.frame(ID=rownames(resDESeq2),
#                        meanA=rowMeans(NormCount[,which(samples_conditions$CONDITIONS == conditionA)]),
#                        meanB=rowMeans(NormCount[,which(samples_conditions$CONDITIONS == conditionB)]),
#                        padj=resDESeq2$padj,
#                        log2FC=resDESeq2$log2FoldChange,
#                        NormCount)                        


#head(diff_counts)

write.table(diff_counts,file=output,
            sep="\t",
            quote=FALSE,
            col.names = TRUE,
            row.names = FALSE)




#samples_conditions[samples_conditions$CONDITIONS == "RC" , "sample"]

#conditions_ordre <- colData[colData == "RC"]

#conditions_ordre
# Sélection des colonnes correspondant à la condition "RC"
#liste=c("RC2","RC3")
#NormCount <- as.matrix(NormCount)
#colonnes_RC <- NormCount[,liste]
#head(colonnes_RC)
#NormCount[,subset(colData, condition == conditionA)]


#meanA=rowMeans(NormCount[,colnames(subset(colData, condition == conditionA))])
#meanA

#head(NormCount$pvalue)
#names(pvalueAll)  = c("tag","pvalue")
#adjPvalue         = p.adjust(as.numeric(as.character(pvalueAll[,"pvalue"])),"BH")
#data=data.frame(ID=rownames(resDESeq2),
#                         meanA=rowMeans(NormCount[,rownames(subset(colData, condition == conditionA))]),
#                         meanB=rowMeans(NormCount[,rownames(subset(colData, condition == conditionB))]),
#                         log2FC=resDESeq2$log2FoldChange,
#                         NormCount)
#
#head(data)
  #COLLECT COUNTS
  
  #names(NormCount) <- NormCount_names
  
  # WRITE A TSV WITH THIS FORMAT FOR THE CURRENT CHUNK
  # Kmer_ID
  # meanA
  # meanB
  # log2FC
  # NormCount
  
#  write.table(data.frame(ID=rownames(resDESeq2),
#                         meanA=rowMeans(NormCount[,rownames(subset(colData, condition == conditionA))]),
#                         meanB=rowMeans(NormCount[,rownames(subset(colData, condition == conditionB))]),
#                         log2FC=resDESeq2$log2FoldChange,
#                         NormCount),
#              file=gzfile(paste(output_tmp_DESeq2,i,"_dataDESeq2_part_tmp.gz", sep="")),
#              sep="\t",quote=FALSE,
#              row.names = FALSE,
#              col.names = FALSE)
#  
#  # WRITE PVALUES FOR THE CURRENT CHUNK
#  write.table(data.frame(ID=rownames(resDESeq2),pvalue=resDESeq2$pvalue),
#              file=gzfile(paste(output_tmp_DESeq2,i,"_pvalue_part_tmp.gz",sep="")),
#              sep="\t",quote=FALSE,
#              row.names = FALSE,
#              col.names = FALSE)
#  
#  # Remove processed chunk
#  system(paste("rm",lst_files[i]))
#  
#}) #END FOREACH
#
#system(paste("rm -rf", output_tmp_chunks))
#
#logging("Foreach done")
#
##MERGE ALL CHUNKS PVALUE INTO A FILE
#system(paste("find", output_tmp_DESeq2, "-name '*_pvalue_part_tmp.gz' | xargs cat >", output_pvalue_all))
#
#logging(paste("Pvalues merged into",output_pvalue_all))
#
##MERGE ALL CHUNKS DESeq2 INTO A FILE
#system(paste("find", output_tmp_DESeq2, "-name '*_dataDESeq2_part_tmp.gz' | xargs cat >", dataDESeq2All))
#
#logging(paste("DESeq2 results merged into",dataDESeq2All))
#
## REMOVE DESeq2 CHUNKS RESULTS
#system(paste("rm -rf", output_tmp_DESeq2))
#
##CREATE AND WRITE THE ADJUSTED PVALUE UNDER THRESHOLD WITH THEIR ID
#pvalueAll         = read.table(output_pvalue_all, header=F, stringsAsFactors=F)
#names(pvalueAll)  = c("tag","pvalue")
#adjPvalue         = p.adjust(as.numeric(as.character(pvalueAll[,"pvalue"])),"BH")
#
#adjPvalue_dataframe = data.frame(tag=pvalueAll$tag,
#                                 pvalue=adjPvalue)
#
#write.table(adjPvalue_dataframe,
#            file=gzfile(adj_pvalue),
#            sep="\t",
#            quote=FALSE,
#            col.names = FALSE,
#            row.names = FALSE)
#
#logging("Pvalues are adjusted")
#
##LEFT JOIN INTO dataDESeq2All
##GET ALL THE INFORMATION (ID,MEAN_A,MEAN_B,LOG2FC,COUNTS) FOR DE KMERS
#system(paste("echo \"LANG=en_EN join <(zcat ", adj_pvalue," | LANG=en_EN sort -n -k1) <(zcat ", dataDESeq2All," | LANG=en_EN sort -n -k1) | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs(\\$5) >=", log2fc_threshold, " && \\$2 <= ", pvalue_threshold, ") print \\$0}' | tr ' ' '\t' | gzip > ", dataDESeq2Filtered, "\" | bash", sep=""))
#system(paste("rm", adj_pvalue, dataDESeq2All))
#
#logging("Get counts for pvalues that passed the filter")
#
##CREATE THE FINAL HEADER USING ADJ_PVALUE AND DATADESeq2ALL ONES AND COMPRESS THE FILE
## CREATE THE HEADER FOR THE DESeq2 TABLE RESULT
#
##SAVE THE HEADER
#system(paste("echo 'tag\tpvalue\tmeanA\tmeanB\tlog2FC' | paste - ", header_kmer_counts," | gzip > ", output_diff_counts))
#system(paste("cat", dataDESeq2Filtered, ">>", output_diff_counts))
#system(paste("rm", dataDESeq2Filtered))
#
#logging("Analysis done")