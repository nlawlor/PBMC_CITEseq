# pathway analysis of doublets removed
rm(list = ls())
library(Matrix)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

mono_lps <- read.csv("../../Monocyte.doublets.removed.LPS.vs.Unactivated.marker.genes.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)

# separate out positive and negative genes into lists
lps_up <- mono_lps[mono_lps$avg_logFC > 0.25 & mono_lps$p_val_adj < 0.05, ]
lps_dw <- mono_lps[mono_lps$avg_logFC < -0.25 & mono_lps$p_val_adj < 0.05, ]

# convert to entrez ids
lps.gene.df <- bitr(rownames(lps_up), fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)
lps.gene.dw.df <- bitr(rownames(lps_dw), fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)

id_lps <- which(lps.gene.df$SYMBOL %in% rownames(lps_up))
id_lps_dw <- which(lps.gene.dw.df$SYMBOL %in% rownames(lps_dw))

mydf <- data.frame(Entrez = c(lps.gene.df$ENTREZID[id_lps], lps.gene.dw.df$ENTREZID[id_lps_dw]),
                   group = c(rep("LPS Induced", length(id_lps)), rep("LPS Reduced", length(id_lps_dw))))

# GO pathways
formula_res <- compareCluster(Entrez~group, data=mydf, fun="enrichGO", OrgDb='org.Hs.eg.db')
print(dotplot(formula_res, showCategory = 20) + ggtitle(paste("GO", sep = " ")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(size = 10)))


