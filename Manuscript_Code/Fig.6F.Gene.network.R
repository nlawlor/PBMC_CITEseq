# pathway analysis of doublets removed
rm(list = ls())
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# load in marker genes
mono_lps <- read.csv("../../Monocyte.doublets.removed.LPS.vs.Unactivated.marker.genes.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)

# separate out positive and negative genes into lists
lps_up <- mono_lps[abs(mono_lps$avg_logFC) > 0.25 & mono_lps$p_val_adj < 0.05, ]

# convert to entrez ids
lps.gene.df <- bitr(rownames(lps_up), fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)

# match up fc values with genes
id_fc <- match(lps.gene.df$SYMBOL, table = rownames(lps_up))
lps.gene.df$FC <- lps_up$avg_logFC[id_fc]
lps.gene.df$FDR <- lps_up$p_val_adj[id_fc]

# make a vector of entrez ids and fc
genes_use <- lps.gene.df$FC
names(genes_use) <- lps.gene.df$ENTREZID

ekegg <- enrichKEGG(names(genes_use), organism = "hsa", qvalueCutoff = 0.05)
## convert gene ID to Symbol
edox <- setReadable(ekegg, 'org.Hs.eg.db', 'ENTREZID')
edox_sel <- edox
# subset the pathway results
edox_sel@result <- edox_sel@result[edox_sel@result$Description %in% c("NOD-like receptor signaling pathway", "Toll-like receptor signaling pathway",
                                                                      "Phagosome", "Graft-versus-host disease", 
                                                                      "NF-kappa B signaling pathway"),]
# plot by p-value of each gene in diff exp analysis
genes_use <- -log10(lps.gene.df$FDR)
id_inf <- which(is.infinite(genes_use))
genes_use[id_inf] <- 300
names(genes_use) <- lps.gene.df$ENTREZID

# identify which genes were reduced and add a signage
id_nge <- which(lps.gene.df$FC < 0)
genes_use[id_nge] <- -genes_use[id_nge]

p1 <- cnetplot(edox_sel, foldChange=genes_use, showCategory = 5, colorEdge = FALSE) + 
  ggtitle("KEGG Pathways Monocyte-LPS Response Genes") + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(color = "-log10 FDR of Gene",
       size = "# of Genes")
p1

