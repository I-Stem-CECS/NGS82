# 2023-03-01 NGS82 
# RNAseq GSDIII
# H. Polveche

library(tidyverse)
library(cowplot)

#######################
### Violin Plot AGL ###
#######################

cts <- read.csv2("../results/GeneralView/NGS82_readscounts_norm_19311genes.csv", 
                 row.names = 1, dec = ".")
coldata <- read.csv2("./coldata.csv", sep = "\t")

AGL <- cts[which(cts$symbol %in% "AGL"),] # ENSG00000162688
AGL <- AGL[, c(2:ncol(AGL))]

AGL <- AGL %>% 
  pivot_longer(everything(), names_to = "samples", values_to = "count") 

AGL.m <- as_tibble(merge(AGL, coldata, by.x = "samples", by.y = "LibQ", all = T))

png("../results/NGS82_ViolinPlot_GSDIIIvsWT_CellLines.png", res = 300, 
    width = 1500, height = 1000)
ggplot(AGL.m, aes(x = cell_line, y = log10(count + 1), fill = condition)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill = "black") +
  theme_cowplot() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         axis.title.y = element_text(size = 10)) +
  labs(fill = "AGL") + 
  ylab("Normalized counts (log10)") +
  xlab("Cell line") +
  scale_fill_manual(values=c("firebrick", "steelblue"))
dev.off()

png("../results/NGS82_ViolinPlot_GSDIIIvsWT.png", res = 300, 
    width = 1000, height = 1000)
ggplot(AGL.m, aes(x = condition, y = log10(count + 1), fill = condition)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill = "black") +
  theme_cowplot() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position = "none",
         axis.title.y = element_text(size = 10)) +
  ylab("Normalized counts (log10)") +
  xlab(" ") +
  scale_fill_manual(values=c("firebrick", "steelblue"))
dev.off()
