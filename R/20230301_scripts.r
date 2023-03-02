# 2023-03-01 NGS82 
# RNAseq GSDIII

library(RColorBrewer)
library(gplots)
library(DESeq2)
library("pheatmap")
library("ggplot2")
library("gridExtra")
library(reshape)

library("tidyverse")
library("ggforce")
library("cowplot")
library("reshape2")
library(viridis)

directory <- "../COUNTSTAR/"

sampleFiles <- grep("*txt",list.files(directory),value=TRUE)
sampleCondition <- sub("*_htseq_count.txt","",sampleFiles)
sampleTable <- data.frame(sampleName = sampleCondition,
                          fileName = sampleFiles,
                          condition = sampleCondition)

write.csv2(sampleTable, "../results/sampleTable.csv")

names_genes <- read.csv2("../../genomes/GRCh37.87/idGene2Symbole_GRCh37.87_ISTEM_WITHOUT_uncharacterised.txt", sep=",", header=F , stringsAsFactors=F)
colnames(names_genes) <- c("id_ensembl","name_gene")

coldata <- read.csv2("./coldata.csv", sep="\t", header=T)
rownames(coldata) <- coldata$sampleName

coldata$condition <- as.factor(coldata$condition)
coldata$cell_line <- as.factor(coldata$cell_line)

sampleTable2 <- merge(sampleTable[,1:2], coldata, by.x="sampleName", by.y="LibQ")
colnames(sampleTable2)[3] <- "Name"
rownames(sampleTable2) <- sampleTable2$sampleName


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2,
                                       directory = directory,
				       design = ~ condition)

dds <- DESeq(ddsHTSeq)

cts <- counts(dds, norm=F) #57905
cts.0 <- cts[which(rownames(cts) %in% names_genes$id_ensembl),] # 19311

counts_raw <- data.frame(symbol=NA, cts.0)
for (i in 1:nrow(counts_raw) ){
  counts_raw[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_raw)[i]),2]
}
write.table(counts_raw,file="../results/NGS82_readscounts_raw_19311genes.csv", sep=";")

cts.0.norm <- counts(dds, norm=T)
cts.0.norm <- cts.0.norm[which(rownames(cts) %in% names_genes$id_ensembl),]

counts_normalise_HP <- data.frame(symbol=NA, cts.0.norm)
for (i in 1:nrow(counts_normalise_HP) ){
  counts_normalise_HP[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_normalise_HP)[i]),2]
}
write.table(counts_normalise_HP,file="../results/NGS82_readscounts_norm_19311genes.csv", sep=";")

cts.1 <- data.frame(cts.0, moy=NA)
for (i in 1:nrow(cts.1)){
  cts.1[i, "moy"] <- mean(as.numeric(cts.1[i , c(1:ncol(cts))]))
}
cts.2 <- cts.1[which(cts.1$moy > 5),]
cts.filtre <- cts.2[, c(1:ncol(cts))] # 12484
colnames(cts.filtre) <- colnames(cts)

##################################
####### dds cts.filtre coldata ###
##################################

for (i in 1:ncol(cts.filtre)){
  colnames(cts.filtre)[i] <- coldata[which(coldata$LibQ %in% colnames(cts.filtre)[i]), "sampleName"]
}

cts.filtre <- cts.filtre[,c(order(colnames(cts.filtre)))]
coldata <- coldata[order(coldata$sampleName),]

dds <- DESeqDataSetFromMatrix(countData = cts.filtre, colData = coldata,
                              design = ~ condition )
dds <- DESeq(dds)

####################################
## Comptages bruts et Normalisés ##
####################################

#### Bruts
counts_raw_HP <- counts(dds, norm=F)

counts_raw <- data.frame(symbol=NA, counts_raw_HP)
for (i in 1:nrow(counts_raw) ){
        counts_raw[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_raw)[i]),2]
}

write.table(counts_raw,file="../results/NGS82_readscounts_raw_12484genes.csv", sep=";")

#### Norm
counts_normalise <- counts(dds, norm=T)
counts_normalise_HP <- data.frame(symbol=NA, counts_normalise)
for (i in 1:nrow(counts_normalise) ){
        counts_normalise_HP[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_normalise_HP)[i]),2]
}

write.table(counts_normalise_HP,file="../results/NGS82_readscounts_norm_12484genes.csv", sep=";", dec=",")


#### GENERAL VIEW ###

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)

pdf("../results/plots_general_view_NGS82.pdf")

condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
condcols <- condcols[1:length(unique(coldata$condition))]
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[as.factor(coldata$condition)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[as.factor(coldata$condition)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')

pcaData <- plotPCA(rld, intgroup=c("condition", "cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, shape=condition, color=cell_line)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15), cexRow=0.5, cexCol=0.5)


plotDispEsts(dds)


dev.off()



###########################
### Counts observations ###
###########################

# counts heatmap


annotation_col <- data.frame( condition = coldata$condition, cell_line = coldata$cell_line) 
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue"), cell_line = brewer.pal(8, "Paired")) 
names(mat_colors$condition) <- unique(coldata$condition)
names(mat_colors$cell_line) <- unique(coldata$cell_line)

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#cts.NGS66.2 <- cts.NGS66[,c(2:ncol(cts.NGS66))]
cts.norm.df <- as.data.frame(counts_normalise)
png("../results/NGS82_pheatmap_counts_12484.png" , width = 1000, height = 1000)
pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, 
          cluster_cols=TRUE, annotation_col = annotation_col ,
          annotation_colors = mat_colors , angle_col = "45")
dev.off()


png("../results/NGS82_pheatmap_counts_12484_BLUE.png" , width = 1000, height = 1000)
pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = colours,
          cluster_cols=TRUE, annotation_col = annotation_col ,
          annotation_colors = mat_colors , angle_col = "45")
dev.off()

png("../results/NGS82_pheatmap_counts_12484_VIRIDIS.png" , width = 1000, height = 1000)
pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col ,
          annotation_colors = mat_colors , angle_col = "45")
dev.off()


#####################
### GSDIII vs WT  ###
#####################

cdition1 <- "GSDIII"
cdition2Ctrl <- "WT"
contrastO <- "condition"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre, 
		colData = coldata , 
		design = ~ condition)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 9
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 3
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 6

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS82.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
condcols <- condcols[1:length(unique(coldata$condition))]
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("condition","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, shape=condition, color=cell_line)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
     )      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
       )
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( condition = coldata$condition, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(8, "Dark2")
                   )
names(mat_colors$condition) <- c(cdition1, cdition2Ctrl)
names(mat_colors$cell_line) <- unique(coldata$cell_line)

pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## 


## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
	mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$condition %in% cdition1) | 
                         (coldata$condition %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$condition))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)

for (i in 1:nrow(finalDE) ){
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/NGS82_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/NGS82_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))

#####
#####
#####

coldata$cell_line <- as.character(coldata$cell_line)
coldata[which(coldata$condition %in% "WT"),"cell_line"]  <- c(rep("WT", 9 ))
coldata$cell_line <- as.factor(coldata$cell_line)

##########################
### GSDIII-0576 vs WT  ###
##########################

cdition1 <- "0576c12"
cdition2Ctrl <- "WT"
contrastO <- "cell_line"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$cell_line %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$cell_line %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$cell_line %in% cdition1) | 
                                                        (coldata$cell_line %in% cdition2Ctrl) ),] , 
                              design = ~ cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 171
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 101
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 70

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS82.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$cell_line)), name = 'Paired')
#condcols <- condcols[1:length(unique(coldata$cell_line))]
names(condcols)=unique(coldata$cell_line)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=cell_line)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
)      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
)
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(
                   cell_line = brewer.pal(length(unique(coldata$cell_line)), "Dark2")
)
names(mat_colors$cell_line) <- unique(coldata$cell_line)

pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## 


## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
  mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$cell_line %in% cdition1) | 
                         (coldata$cell_line %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$cell_line))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)

for (i in 1:nrow(finalDE) ){
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/NGS82_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/NGS82_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))

