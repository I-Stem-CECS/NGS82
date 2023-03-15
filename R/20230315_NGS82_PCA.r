##############################
### 2023-03-09 PCA NGS82  ####
### H. Polveche ##############
##############################

# Source :
# http://www.sthda.com/french/articles/38-methodes-des-composantes-principales-dans-r-guide-pratique/73-acp-analyse-en-composantes-principales-avec-r-l-essentiel/

library("FactoMineR") # PCA
library("factoextra") # PCA
library("ggplot2")    # Graphiques
library(RColorBrewer) # Couleurs

cts <- read.csv2("../results/GeneralView/NGS82_readscounts_norm_19311genes.csv", 
                 row.names = 1, dec = ".")


rownames(cts) <- cts$symbol
cts <- cts[,c(2:ncol(cts))]

coldata <- read.csv2("./coldata.csv", sep = "\t")
colnames(cts) <- coldata$sampleName

res.pca <- PCA(t(cts), graph=F, scale.unit = T )
eig.val <- get_eigenvalue(res.pca)

png("../results/PCA/NGS82_eig_val.png", res = 300, width = 1000, height = 1000)
  fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 20))
dev.off()

write.csv2(eig.val, "../results/PCA/NGS82_eig_val.csv")

var <- get_pca_var(res.pca)
var

summary(var$contrib[,"Dim.1"])
fviz_contrib(res.pca, choice = "var", axes = 1, top = 20)

png("../results/PCA/NGS82_AGL_contrib.png", res = 300, width = 2000, height = 2000)
fviz_pca(res.pca,  select.var = list(name = c("AGL")),
         col.ind = "cos2", axes = c(3, 4), gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) +
  labs(title ="genes")
dev.off()

res.desc <- dimdesc(res.pca, axes = c(1:5), proba = 0.05)
str(res.desc)

write.csv2(na.omit(res.desc$Dim.1$quanti), "../results/PCA/NGS82_Dim1_Genes_contrib.csv")
write.csv2(na.omit(res.desc$Dim.2$quanti), "../results/PCA/NGS82_Dim2_Genes_contrib.csv")
write.csv2(na.omit(res.desc$Dim.3$quanti), "../results/PCA/NGS82_Dim3_Genes_contrib.csv")
write.csv2(na.omit(res.desc$Dim.4$quanti), "../results/PCA/NGS82_Dim4_Genes_contrib.csv")
write.csv2(na.omit(res.desc$Dim.5$quanti), "../results/PCA/NGS82_Dim5_Genes_contrib.csv")

ind <- get_pca_ind(res.pca)
ind

ind.coord <- ind$coord
write.csv2(ind.coord, "../results/PCA/NGS82_ind_coord.csv")


png("../results/PCA/NGS82_pca_ind_PC1-2.png", res = 300, width = 2000, height = 2000)
fviz_pca_ind(res.pca, col.ind = "cos2",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE # Évite le chevauchement de texte
)
dev.off()


png("../results/PCA/NGS82_pca_ind_PC3-4.png", res = 300, width = 2000, height = 2000)
fviz_pca_ind(res.pca, col.ind = "cos2", axes = c(3,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
)
dev.off()

# Contribution totale sur PC1 et PC2
png("../results/PCA/NGS82_contrib_PC1.png", res = 300, width = 2000, height = 2000)
fviz_contrib(res.pca, choice = "ind", axes = 1)
dev.off()
# Contribution totale sur PC3 et PC4
png("../results/PCA/NGS82_contrib_PC2.png", res = 300, width = 2000, height = 2000)
fviz_contrib(res.pca, choice = "ind", axes = 2)
dev.off()

png("../results/PCA/NGS82_contrib_PC3.png", res = 300, width = 2000, height = 2000)
fviz_contrib(res.pca, choice = "ind", axes = 3)
dev.off()

png("../results/PCA/NGS82_contrib_PC4.png", res = 300, width = 2000, height = 2000)
fviz_contrib(res.pca, choice = "ind", axes = 4)
dev.off()

## Dim.3 - correl avec AGL ?

rownames(coldata) <- coldata$sampleName
ind.coldata <- data.frame(coldata,ind.coord)

png("../results/PCA/NGS82_cPCA_1-2.png", res = 300, width = 2000, height = 2000)
p <- ggplot(data=ind.coldata, aes(x=Dim.1, y=Dim.2, colour=cell_line, shape = condition))
p <- p + geom_point(size=4)
print(p)
dev.off()

png("../results/PCA/NGS82_cPCA_3-4.png", res = 300, width = 2000, height = 2000)
p2 <- ggplot(data=ind.coldata, aes(x=Dim.3, y=Dim.4, colour=cell_line, shape = condition))
p2 <- p2 + geom_point(size=4)
print(p2)
dev.off()
