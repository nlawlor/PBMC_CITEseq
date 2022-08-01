library(Matrix)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(cowplot)

# load in metadata with scores
meta <- read.csv("../../Meta_CD4_CD8_NK_05122020.csv", header = T, check.names = F, stringsAsFactors = F,
                 row.names = 1)

# read in seurat object
ser <- readRDS("../../CD8T.CD4T.and.NK.cells.no.doublets.seurat.Rds")
rownames(meta) <- gsub(x = rownames(meta), pattern = "-0", replacement = "")
rownames(meta) <- gsub(x = rownames(meta), pattern = "-1", replacement = "")
rownames(meta) <- gsub(x = rownames(meta), pattern = "-2", replacement = "")
table(colnames(ser) %in% rownames(meta))
table(colnames(ser) == rownames(meta))

# order metadata
meta_ord <- meta[colnames(ser), ]
table(colnames(ser) == rownames(meta_ord))

# add zeros to mising
id_1 <- which(is.na(meta_ord$NK_scores))
id_2 <- which(is.na(meta_ord$Temra_scores))
meta_ord$Temra_scores[id_2] <- 0
id_3 <- which(is.na(meta_ord$MAIT_scores))
meta_ord$MAIT_scores[id_3] <- 0
id_4 <- which(is.na(meta_ord$Senes_scores))


ser <- AddMetaData(ser, metadata = meta_ord$MAIT_scores, col.name = "MAIT_Score")
ser <- AddMetaData(ser, metadata = meta_ord$NK_scores, col.name = "NK_Score")
ser <- AddMetaData(ser, metadata = meta_ord$Temra_scores, col.name = "TEMRA_Score")
ser <- AddMetaData(ser, metadata = meta_ord$Senes_scores, col.name = "Senes_Score")

# reorder the data by cell types and subsets
levels(ser$Cell_Type)
ser$Cell_Type <- factor(ser$Cell_Type, levels = c("Actv_NK", "NK",
                                                  "CD4T_Naive", "CD4T_Naive_Actv", "CD4T_Mem", "CD4T_Mem_Actv", 
                                                  "CCR5+CCR7+", "CXCR5+", "Th17",
                                                  "CD8T_Naive", "CD8T_Naive_Actv", "CD8T_Mem", "CD8T_Mem_Actv", 
                                                  "CD57+", "MAIT"))

dp1 <- VlnPlot(ser, features = "MAIT_Score", group.by = "Cell_Type",
               pt.size = 0) + NoLegend() + theme(axis.text.x = element_blank())
dp2 <- VlnPlot(ser, features = "NK_Score", group.by = "Cell_Type",
               pt.size = 0) + NoLegend() + theme(axis.text.x = element_blank())
dp3 <- VlnPlot(ser, features = "Senes_Score", group.by = "Cell_Type",
               pt.size = 0) + NoLegend()

cowplot::plot_grid(dp1, dp2, dp3, ncol = 1, align = "v")

dp3 <- VlnPlot(ser, features = "Senes_Score", group.by = "Cell_Type",
               pt.size = 0) + NoLegend() + theme(axis.text.x = element_blank())
cowplot::plot_grid(dp1, dp2, dp3, ncol = 1, align = "v")

