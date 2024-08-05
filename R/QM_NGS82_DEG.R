library("dplyr")
library("DESeq2")
library("tidyverse")
library("biomaRt")
library("openxlsx")

folder_exp <- "data/"

data <- read.csv2(paste0(folder_exp, "NGS82_readscounts_raw_12484genes.csv")) %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::select(-symbol)

transform_ensembl <- function(dataset) {
  ensembl <- useEnsembl(biomart = "genes")
  ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
  
  test <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
    mart = ensembl
  )
  
  dataset <- left_join(data, test, by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
    filter(gene_biotype == "protein_coding") %>%
    filter(hgnc_symbol != "") %>%
    dplyr::select(-gene_biotype, -Ensembl_ID) %>%
    dplyr::rename(symbol = hgnc_symbol)
  
  return(dataset)
}

data <- transform_ensembl(data)

data <- data %>% column_to_rownames(var = "symbol")


coldata <- read.xlsx(xlsxFile = paste0(folder_exp, "resume_NGS82.xlsx"), sheet = 3) %>%
  dplyr::select(sampleName, condition, cell_line) %>%
  mutate(sampleName = str_remove(sampleName, "MyoT_starv2D_")) %>%
  dplyr::filter(cell_line != "AGL_KO") %>%
  mutate_all(as.factor)

data <- data %>% mutate_all(as.numeric)

dds <- DESeqDataSetFromMatrix(
  countData = data[, as.character(coldata$sampleName)],
  colData = coldata,
  design = ~condition
)

dds$condition <- relevel(dds$condition, ref = "WT")

keep <- rowMeans(counts(dds)) >= 5
dds <- dds[keep,]
dds <- DESeq(dds)

res <- results(dds, altHypothesis = "greaterAbs", alpha = 0.05) %>%
  subset(baseMean > 40)

res_df <- res %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > log2(2)) %>%
  dplyr::rename(symbol = ID)

