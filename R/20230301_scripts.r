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

directory <- "../COUNTSTAR//"

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
cts.filtre <- cts.2[, c(1:ncol(cts))] # 10344
colnames(cts.filtre) <- colnames(cts)

##################################
####### dds cts.filtre coldata ###
##################################

cts.filtre <- cts.filtre[,c(order(colnames(cts.filtre)))]
coldata <- coldata[order(coldata$sampleName),]

dds <- DESeqDataSetFromMatrix(countData = cts.filtre, colData = coldata,
                              design = ~ state + cell_line )
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

write.table(counts_raw,file="./results/NGS78_readscounts_raw_10311genes.csv", sep=";")

#### Norm
counts_normalise <- counts(dds, norm=T)
counts_normalise_HP <- data.frame(symbol=NA, counts_normalise)
for (i in 1:nrow(counts_normalise) ){
        counts_normalise_HP[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_normalise_HP)[i]),2]
}

write.table(counts_normalise_HP,file="./results/NGS78_readscounts_norm_10311genes.csv", sep=";", dec=",")


#### GENERAL VIEW ###

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)

pdf("./results/plots_general_view_NGS78.pdf")

condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[as.factor(coldata$state)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[as.factor(coldata$state)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')

pcaData <- plotPCA(rld, intgroup=c("state", "cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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


##################################
####### dds without d385_P0 et d384_J6 d383_J13 ###
##################################

cts.filtre <- cts.filtre[,c(order(colnames(cts.filtre)))]
coldata <- coldata[order(coldata$sampleName),]

cts.filtre <- cts.filtre[,c(1:16, 18:ncol(cts.filtre))]
coldata <- coldata[c(1:16, 18:nrow(coldata)),]

cts.filtre <- cts.filtre[,c(1:9, 11:ncol(cts.filtre))]
coldata <- coldata[c(1:9, 11:nrow(coldata)),]

cts.filtre <- cts.filtre[,c(1, 3:ncol(cts.filtre))]
coldata <- coldata[c(1, 3:nrow(coldata)),]


dds <- DESeqDataSetFromMatrix(countData = cts.filtre, colData = coldata,
                              design = ~ state + cell_line )
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

write.table(counts_raw,file="./results/NGS78_readscounts_raw_10311genes.csv", sep=";")

#### Norm
counts_normalise <- counts(dds, norm=T)
counts_normalise_HP <- data.frame(symbol=NA, counts_normalise)
for (i in 1:nrow(counts_normalise) ){
  counts_normalise_HP[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_normalise_HP)[i]),2]
}

write.table(counts_normalise_HP,file="./results/NGS78_readscounts_norm_10311genes.csv", sep=";", dec=",")


#### GENERAL VIEW ###

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)

pdf("./results/plots_general_view_NGS78.pdf")

condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[as.factor(coldata$state)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[as.factor(coldata$state)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')

pcaData <- plotPCA(rld, intgroup=c("state", "cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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


annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line) #, genotype = coldata$genotype )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = brewer.pal(6, "Spectral"), cell_line = brewer.pal(4, "Paired")) #  , genotype = c("mediumpurple", "hotpink4"))
names(mat_colors$state) <- unique(coldata$state)
names(mat_colors$cell_line) <- unique(coldata$cell_line)

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#cts.NGS66.2 <- cts.NGS66[,c(2:ncol(cts.NGS66))]
cts.norm.df <- as.data.frame(counts_normalise)
png("./results/NGS78_pheatmap_counts_10344.png" , width = 1000, height = 1000)
pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, 
          cluster_cols=TRUE, annotation_col = annotation_col ,
          annotation_colors = mat_colors , angle_col = "45")
dev.off()


png("./results/NGS78_pheatmap_counts_10344_BLUE.png" , width = 1000, height = 1000)
pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = colours,
          cluster_cols=TRUE, annotation_col = annotation_col ,
          annotation_colors = mat_colors , angle_col = "45")
dev.off()




#################
### D6 vs D1  ###
#################

cdition1 <- "D6"
cdition2Ctrl <- "D1"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
              which((colnames(cts.filtre) %in% coldata[ which(
                coldata$state %in% cdition1), "sampleName"]) | 
		(colnames(cts.filtre) %in% coldata[ which(
		  coldata$state %in% cdition2Ctrl ), "sampleName"])
		)], 
		colData = coldata[which((coldata$state %in% cdition1) | 
		                          (coldata$state %in% cdition2Ctrl) ),] , 
		design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 14
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 14
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
                   )
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-j1_S7`=NA, `d384-j1_S13`=NA, `d385-j1_S19`=NA, `dk7-j1_S1`=NA, 
                      `d383-j6_S8`=NA, `d385-j6_S20`=NA, `dk7-j6_S2`=NA
                      )

for (i in 1:nrow(finalDE) ){
	##D1
	finalDE[i,c(7:10)] <- counts_norm_HP[which(
	                rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
	                c(colnames(counts_norm_HP[, c(1,3,4,6)]))]
	##D6
	finalDE[i,c(11:13)] <- counts_norm_HP[which(
	                rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
	                c(colnames(counts_norm_HP[, -c(1,3,4,6)]))]

	if (is.na(finalDE[i,"log2FoldChange"])) {
		finalDE[i, "differential_expression"] <- NA
	} else if (finalDE[i, "log2FoldChange"] >= 0 ) {
		finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
	} else {
		finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
	}
	
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))




#################
### D13 vs D1  ###
#################

cdition1 <- "D13"
cdition2Ctrl <- "D1"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 22
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 22
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-j1_S7`=NA, `d384-j1_S13`=NA, `d385-j1_S19`=NA, `dk7-j1_S1`=NA, 
                      `d384-j13_S15`=NA, `d385-j13_S21`=NA, `dk7-j13_S3`=NA
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:10)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(1,2,4,6)]))]
  ##D6
  finalDE[i,c(11:13)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(1,2,4,6)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))



#################
### D17 vs D1  ###
#################

cdition1 <- "D17"
cdition2Ctrl <- "D1"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 22
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 22
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

head(counts_norm_HP)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-j1_S7`=NA, `d384-j1_S13`=NA, `d385-j1_S19`=NA, `dk7-j1_S1`=NA, 
                      `d383-j17_S10`=NA, `d384-j17_S16`=NA, `d385-j17_S22`=NA, `dk7-j17_S4`=NA
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:10)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(1,3,5,7)]))]
  ##D6
  finalDE[i,c(11:14)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(1,3,5,7)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))


#################
### P0 vs D1  ###
#################

cdition1 <- "P0"
cdition2Ctrl <- "D1"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 22
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 22
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

head(counts_norm_HP)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-j1_S7`=NA, `d384-j1_S13`=NA, `d385-j1_S19`=NA, `dk7-j1_S1`=NA, 
                      `d383-p0_S11`=NA, `d384-p0_S17`=NA, `dk7-p0_S5`=NA
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:10)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(1,3,5,6)]))]
  ##D6
  finalDE[i,c(11:13)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(1,3,5,6)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))


#################
### P1 vs D1  ###
#################

cdition1 <- "P1"
cdition2Ctrl <- "D1"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 14
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 14
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

head(counts_norm_HP)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-j1_S7`=NA, `d384-j1_S13`=NA, `d385-j1_S19`=NA, `dk7-j1_S1`=NA, 
                      `d383-p1_S12`=NA, `d384-p1_S18`=NA, `d385-p1_S24`=NA, `dk7-p1_S6`=NA
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:10)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(1,3,5,7)]))]
  ##D6
  finalDE[i,c(11:14)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(1,3,5,7)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))




#################
### D13 vs D6  ###
#################

cdition1 <- "D13"
cdition2Ctrl <- "D6"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 1
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 1

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

head(counts_norm_HP)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-j6_S8`=NA, `d385-j6_S20`=NA, `dk7-j6_S2`=NA,
                      `d384-j13_S15`=NA, `d385-j13_S21`=NA, `dk7-j13_S3`=NA
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:9)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(1,4,6)]))]
  ##D6
  finalDE[i,c(10:12)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(1,4,6)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))




#################
### D17 vs D13  ###
#################

cdition1 <- "D17"
cdition2Ctrl <- "D13"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

head(counts_norm_HP)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d384-j13_S15`=NA, `d385-j13_S21`=NA, `dk7-j13_S3`=NA,
                      `d383-j17_S10`=NA, `d384-j17_S16`=NA, `d385-j17_S22`=NA, `dk7-j17_S4`=NA
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:9)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(2,4,6)]))]
  ##D6
  finalDE[i,c(10:13)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(2,4,6)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))




#################
### P0 vs D17 ###
#################

cdition1 <- "P0"
cdition2Ctrl <- "D17"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 3
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 3
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

head(counts_norm_HP)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-j17_S10`=NA, `d384-j17_S16`=NA, `d385-j17_S22`=NA, `dk7-j17_S4`=NA,
                      `d383-p0_S11`=NA, `d384-p0_S17`=NA, `dk7-p0_S5`=NA
                      
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:10)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(1,3,5,6)]))]
  ##D6
  finalDE[i,c(11:13)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(1,3,5,6)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))


#################
### P1 vs P0  ###
#################

cdition1 <- "P1"
cdition2Ctrl <- "P0"
contrastO <- "state"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$state %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$state %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$state %in% cdition1) | 
                                                        (coldata$state %in% cdition2Ctrl) ),] , 
                              design = ~ state + cell_line)



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("./results/plots_",cdition1,"-vs-",cdition2Ctrl,"_NGS78.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$state)), name = 'Paired')
names(condcols)=unique(coldata$state)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("state","cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=state, shape=cell_line)) +
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

annotation_col <- data.frame( state = coldata$state, cell_line = coldata$cell_line )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(state = c("firebrick", "steelblue")  , 
                   cell_line = brewer.pal(4, "Dark2")
)
names(mat_colors$state) <- c(cdition1, cdition2Ctrl)
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

write.csv2(mat_sig, paste0("./results/",cdition1,"-vs-", cdition2Ctrl,"_log2FC0-4.csv"))

## Formatage I-Stem


counts_norm_HP <- counts(dds, norm=T)

head(counts_norm_HP)

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA, 
                      `d383-p0_S11`=NA, `d384-p0_S17`=NA, `dk7-p0_S5`=NA,
                      `d383-p1_S12`=NA, `d384-p1_S18`=NA, `d385-p1_S24`=NA, `dk7-p1_S6`=NA
)

for (i in 1:nrow(finalDE) ){
  ##D1
  finalDE[i,c(7:9)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, c(1,3,6)]))]
  ##D6
  finalDE[i,c(10:13)] <- counts_norm_HP[which(
    rownames(counts_norm_HP) %in% rownames(finalDE)[i]),
    c(colnames(counts_norm_HP[, -c(1,3,6)]))]
  
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

write.csv2(finalDE, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("./results/NGS78_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))



#############################################
### Nombre de genes > 3 reads par samples ###
#############################################

cts.filtre <- cts.filtre[,c(order(colnames(cts.filtre)))]
coldata <- coldata[order(coldata$sampleName),]

dds <- DESeqDataSetFromMatrix(countData = cts.filtre, colData = coldata,
                              design = ~ state + cell_line )
dds <- DESeq(dds)


#### Norm
cts.norm <- counts(dds, norm=T)
cts.TRUE <- as.data.frame(matrix(nc=ncol(cts.norm), nr= nrow(cts.norm)))
rownames(cts.TRUE) <- rownames(cts.norm)
colnames(cts.TRUE) <- colnames(cts.norm)

for (i in 1:nrow(cts.norm)){
  for (j in 1:ncol(cts.norm)){
    if(cts.norm[i, j] > 3){
      cts.TRUE[i, j] <- TRUE
    } else {
      cts.TRUE[i, j] <- FALSE
    }
      
  }
}

aaa <- summary(cts.TRUE)
genes.per.sampes <- t(aaa[2:3,])


write.table(genes.per.sampes,file="./results/genes-per-sampes_supto3reads.csv", sep=";", dec=",")


cts.condition <- as.data.frame(matrix(nc=6, nr = nrow(cts.norm)))
rownames(cts.condition) <- rownames(cts.norm)
colnames(cts.condition) <- c("D1", "D6", "D13","D17", "P0", "P1")

#coldata[which(coldata$state %in% "D1"), "sampleName"]

for (i in 1:ncol(cts.condition)){
  cts.tmp <- cts.norm[, c(coldata[which(coldata$state %in% colnames(cts.condition)[i]), "sampleName"])]
  for (j in 1:nrow(cts.tmp)){
    cts.condition[j, i] <- mean(cts.tmp[j,])
  }
}


cts.condition2 <- as.data.frame(matrix(nc=6, nr = nrow(cts.norm)))
rownames(cts.condition2) <- rownames(cts.norm)
colnames(cts.condition2) <- c("D1", "D6", "D13","D17", "P0", "P1")

for (i in 1:ncol(cts.condition2)){
  print(c(coldata[which(coldata$state %in% colnames(cts.condition2)[i]), "sampleName"]))
  cts.tmp <- cts.TRUE[, c(coldata[which(coldata$state %in% colnames(cts.condition2)[i]), "sampleName"])]
  for (j in 1:nrow(cts.tmp)){
    vec <- as.vector(t(cts.tmp[j,]))
    cts.condition2[j, i] <- table(vec)["TRUE"]
  }
}


cts.condition3 <- data.frame(symbol=NA, cts.condition2)
for (i in 1:nrow(cts.condition3) ){
  cts.condition3[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(cts.condition2)[i]),2]
}
write.table(cts.condition3,file="./results/genesSupTo3_dans_combien_de_replicats.csv", sep=";", dec=",")

for ( i in 1:ncol(cts.condition2)){
  message(colnames(cts.condition2)[i])
  print(table(cts.condition2[,i]))
}


###########################################
#### 2022-05-13 : filtre p-value < 0,05 ###
###########################################

D6 <- read.csv2("./results/D6vsD1/NGS78_D6-vs-D1_log2FC0-4.csv", row.names = 1)
D13 <- read.csv2("./results/D13vsD1/NGS78_D13-vs-D1_log2FC0-4.csv", row.names = 1)
D17 <- read.csv2("./results/D17vsD1/NGS78_D17-vs-D1_log2FC0-4.csv", row.names = 1)
P0 <- read.csv2("./results/P0vsD1/NGS78_P0-vs-D1_log2FC0-4.csv", row.names = 1)
P1 <- read.csv2("./results/P1vsD1/NGS78_P1-vs-D1_log2FC0-4.csv", row.names = 1)

D6.pval <- D6 %>% 
  filter(pvalue <= 0.05 & baseMean >= 20 & abs(log2FoldChange) >= 0.4 )
D13.pval <- D13 %>% 
  filter(pvalue <= 0.05 & baseMean >= 20 & abs(log2FoldChange) >= 0.4 )
D17.pval <- D17 %>% 
  filter(pvalue <= 0.05 & baseMean >= 20 & abs(log2FoldChange) >= 0.4 )
P0.pval <- P0 %>% 
  filter(pvalue <= 0.05 & baseMean >= 20 & abs(log2FoldChange) >= 0.4 )
P1.pval <- P1 %>% 
  filter(pvalue <= 0.05 & baseMean >= 20 & abs(log2FoldChange) >= 0.4 )

write.csv2(D6.pval, 
           "./results/DxvsD1_pval0-05_log2FC_0-4_BM_20/NGS78_D6vsD1_protein_coding_pval0-05_log2FC_0-4_BM_20.csv")

write.csv2(D13.pval, 
           "./results/DxvsD1_pval0-05_log2FC_0-4_BM_20/NGS78_D13vsD1_protein_coding_pval0-05_log2FC_0-4_BM_20.csv")

write.csv2(D17.pval, 
           "./results/DxvsD1_pval0-05_log2FC_0-4_BM_20/NGS78_D17vsD1_protein_coding_pval0-05_log2FC_0-4_BM_20.csv")

write.csv2(P0.pval, 
           "./results/DxvsD1_pval0-05_log2FC_0-4_BM_20/NGS78_P0vsD1_protein_coding_pval0-05_log2FC_0-4_BM_20.csv")

write.csv2(P1.pval, 
           "./results/DxvsD1_pval0-05_log2FC_0-4_BM_20/NGS78_P1vsD1_protein_coding_pval0-05_log2FC_0-4_BM_20.csv")
