###############################
### 2023-03-09 GSEA NGS82  ####
### H. Polveche ###############
###############################
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/


# BiocManager::install("clusterProfiler") #, version = "3.8")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
#http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
library(org.Hs.eg.db)
library(ggnewscale)
library(pathview)
library(DOSE)
library(plyr)
library(R.utils)

  
formatData <- function(path_file){
  file <- as_tibble(read.csv2(path_file, header = T, sep = ";")) %>% 
    filter(padj <= 0.05 & baseMean >= 20)
  file.sig <- file[,c("symbol", "log2FoldChange")]
  list.genes <- file.sig$log2FoldChange
  names(list.genes) <- file.sig$symbol
  # omit any NA values , sort the list in decreasing order (required for clusterProfiler)
  list.genes <- sort(na.omit(list.genes), decreasing = TRUE)
  return(list.genes)
}


#keytypes(org.Hs.eg.db)

gseafun <- function(list.genes, output_path = "./", ont = "BP", nPerm = 10000, minGSSize = 10, 
                    maxGSSize = 800, pvalueCutoff = 1, OrgDb = "org.Hs.eg.db", pAdjustMethod = "BH"){
  message("gseGO")
  gse <- gseGO(geneList= list.genes, 
               ont = ont, 
               keyType = "SYMBOL", 
               #nPerm = nPerm, 
               minGSSize = minGSSize, 
               maxGSSize = maxGSSize, 
               pvalueCutoff = pvalueCutoff, 
               verbose = TRUE, 
               OrgDb = OrgDb, 
               pAdjustMethod = pAdjustMethod)
  
  ego <- enrichGO(gene         = names(list.genes),
                   OrgDb         = OrgDb,
                   keyType       = 'SYMBOL',
                   ont           = ont,
                   pAdjustMethod = pAdjustMethod,
                   pvalueCutoff  = pvalueCutoff,
                   qvalueCutoff  = 1)
  ego <- gofilter(ego, level = 3)
    write.csv2(as.data.frame(gse), paste0(output_path, "gseGO_results.csv"), row.names = F)
    write.csv2(as.data.frame(ego), paste0(output_path, "enriGO_results.csv"), row.names = F)
  
    png(paste0(output_path, "dotplot_enrichGO.png"), width = 2500, height = 3000, res = 300)
    print(dotplot(ego))
    dev.off()
    
    png(paste0(output_path, "cnetplot_enrichGO_padj.png"), width = 2500, height = 2500, res = 300)
    print(cnetplot(ego, categorySize="p.adjust", showCategory = 3))
    dev.off()
    
    gse.sig <- gse[which(gse$p.adjust <= 0.05),]
    
    if (nrow(gse.sig) != 0){   
      message("p-adjust signifcative values")
      message("dotplot")
      # require(DOSE)
      png(paste0(output_path, "dotplot_gseGO_padj.png"), width = 2500, height = 2500, res = 300)
      print(dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign))
      dev.off()
      
      d <- GOSemSim::godata(OrgDb, ont = ont)   
      compare_gse <- enrichplot::pairwise_termsim(gse, semData = d,  method="Wang")
      message("emapplot")
      png(paste0("emapplot_gseGO_padj.png"), width = 2500, height = 2500, res = 300)
      print(emapplot(compare_gse, showCategory = 20))
      dev.off()
      
      # categorySize can be either 'pvalue' or 'geneNum'
      message("gcnetplot")
      png(paste0(output_path, "cnetplot_gseGO_padj.png"), width = 2500, height = 2500, res = 300)
      print(cnetplot(gse, categorySize="p.adjust", foldChange=list.genes, showCategory = 3))
      dev.off()
      
      message("Done.")
      
  } else if (nrow(gse) != 0 & nrow(gse.sig) == 0) { 
      message("no significative (p-adj) values")
      message("dotplot")
      # require(DOSE)
      png(paste0(output_path, "dotplot_gseGO_pvalue.png"), width = 2500, height = 3000, res = 300)
      print(dotplot(gse, showCategory=10, split=".sign", color = "pvalue") + facet_grid(.~.sign))
      dev.off()
      
      d <- GOSemSim::godata(OrgDb, ont = ont)   
      compare_gse <- enrichplot::pairwise_termsim(gse, semData = d,  method="Wang")
      message("emapplot")
      png(paste0(output_path, "emapplot_gseGO_pvalue.png"), width = 2500, height = 2500, res = 300)
      print(emapplot(compare_gse, showCategory = 10, color = "pvalue"))
      dev.off()
      
      # categorySize can be either 'pvalue' or 'geneNum'
      message("gcnetplot")
      png(paste0(output_path, "cnetplot_gseGO_pvalue.png"), width = 2500, height = 2500, res = 300)
      print(cnetplot(gse, categorySize="pvalue", foldChange=list.genes, showCategory = 3))
      dev.off()
      
      message("Done.")
      
    }else {
      message("no enrichment...")
    }
    return(gse)
}


keggfun <- function(liste.genes, output_path = "./", pvalueCutoff = 1, minGSSize = 5,
                    maxGSSize = 800, pAdjustMethod = "BH", nPerm = 10000){

  ENTREZID <-  ldply(mapIds(org.Hs.eg.db, names(liste.genes), 'ENTREZID', 'SYMBOL'))
  liste.entrezid <- liste.genes
  names(liste.entrezid) <- ENTREZID$V1
  
  # Reading KEGG database
  R.utils::setOption("clusterProfiler.download.method","auto")
  message("gseKEGG")
  kk2 <- gseKEGG(geneList     = liste.entrezid,
                 organism     = "hsa",
                 nPerm        = nPerm,
                 minGSSize    = minGSSize,
                 maxGSSize    = maxGSSize,
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = pAdjustMethod,
                 keyType       = "ncbi-geneid")

  
  write.csv2(as.data.frame(kk2), paste0(output_path, "gseKEGG.csv"), row.names = F)
  
  if (nrow(kk2) > 0 ){
    message("dotplot")
      png(paste0(output_path, "dotplot_gseKEGG.png"), width = 2500, height = 3000, res = 300)
      print(dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + 
              facet_grid(.~.sign))
      dev.off()
      
      # # categorySize can be either 'pvalue' or 'geneNum'
      # png(paste0(output_path, "cnetplot_gseKEGG.png"), width = 2500, height = 2500, res = 300)
      # cnetplot(kk2, categorySize="p.adjust", foldChange=list., showCategory = 3)
      # dev.off()
    message("ridgeplot")  
      png(paste0(output_path, "ridgeplot_gseKEGG.png"), width = 2500, height = 3500, res = 300)
      print(ridgeplot(kk2) + labs(x = "enrichment distribution"))
      dev.off()
    
      ### pathview
      dat.kegg <- as.data.frame(kk2) 
      dir.create(paste0(output_path, "KEGG_pathview"))
        # Produce the native KEGG plot (PNG)
      original_path <- getwd()
      setwd(paste0(output_path, "KEGG_pathview/"))
      message("pathview")
      for (i in 1:nrow(dat.kegg)){
        message(i)
        pathway_name <- dat.kegg[i, "ID"]
        dme <- pathview(gene.data=liste.entrezid, 
                        pathway.id=pathway_name, 
                        species = "hsa")  
      }  
      setwd(original_path)
  } else {
      message("no KEGG pathways")
    }
  
    return(kk2)
  }


enrichment_analysis <- function(path_file, output_path){
  genes.sig <- formatData(path_file = path_file)
  gse <- gseafun(genes.sig, output_path = output_path)
  kk <- keggfun(genes.sig, output_path = output_path) 
  return(list("GSEA"=gse, "KEGG"=kk, "list.genes" = genes.sig))
}

path_file <- "../results/GSDIII-vs-WT/NGS82_GSDIII-vs-WT_log2FC0-4_padj0-5_BM20.csv"

genes.sig <- formatData(path_file = path_file)
gseBP <- gseafun(genes.sig, output_path = "../results/GSEA/GSDIII-vs-WT/BP/", ont = "BP")
gseCC <- gseafun(genes.sig, output_path = "../results/GSEA/GSDIII-vs-WT/CC/", ont = "CC")
gseMF <- gseafun(genes.sig, output_path = "../results/GSEA/GSDIII-vs-WT/MF/", ont = "MF")


path_file <- "../results/GSDIII-0303c3-vs-WT/NGS82_0303c3-vs-WT_log2FC0-4_padj0-5_BM20.csv"

genes.sig <- formatData(path_file = path_file)
gseBP <- gseafun(genes.sig, output_path = "../results/GSEA/0303c3-vs-WT/BP/", ont = "BP")
gseCC <- gseafun(genes.sig, output_path = "../results/GSEA/0303c3-vs-WT/CC/", ont = "CC")
gseMF <- gseafun(genes.sig, output_path = "../results/GSEA/0303c3-vs-WT/MF/", ont = "MF")


path_file <- "../results/GSDIII-0576c12-vs-WT/NGS82_0576c12-vs-WT_log2FC0-4_padj0-5_BM20.csv"

genes.sig <- formatData(path_file = path_file)
gseBP <- gseafun(genes.sig, output_path = "../results/GSEA/0576c12-vs-WT/BP/", ont = "BP")
gseCC <- gseafun(genes.sig, output_path = "../results/GSEA/0576c12-vs-WT/CC/", ont = "CC")
gseMF <- gseafun(genes.sig, output_path = "../results/GSEA/0576c12-vs-WT/MF/", ont = "MF")

path_file <- "../results/GSDIII-2523c3-vs-WT/NGS82_2523c3-vs-WT_log2FC0-4_padj0-5_BM20.csv"

genes.sig <- formatData(path_file = path_file)
gseBP <- gseafun(genes.sig, output_path = "../results/GSEA/2523c3-vs-WT/BP/", ont = "BP")
gseCC <- gseafun(genes.sig, output_path = "../results/GSEA/2523c3-vs-WT/CC/", ont = "CC")
gseMF <- gseafun(genes.sig, output_path = "../results/GSEA/2523c3-vs-WT/MF/", ont = "MF")


path_file <- "../results/GSDIII-3390c15-vs-WT/NGS82_3390c15-vs-WT_log2FC0-4_padj0-5_BM20.csv"

genes.sig <- formatData(path_file = path_file)
gseBP <- gseafun(genes.sig, output_path = "../results/GSEA/3390c15-vs-WT/BP/", ont = "BP")
gseCC <- gseafun(genes.sig, output_path = "../results/GSEA/3390c15-vs-WT/CC/", ont = "CC")
gseMF <- gseafun(genes.sig, output_path = "../results/GSEA/3390c15-vs-WT/MF/", ont = "MF")














