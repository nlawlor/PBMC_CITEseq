library(Seurat)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(viridis)

# load in processed seurat object
pbmc <- readRDS("../../CD4T.doublets.removed.no.IgM.seurat.Rds")

# plot by treatment (Fig. 4C)
DimPlot(pbmc, group.by = "Treatment", label = F, cols = c("#984ea3", "#377eb8", "#4daf4a"))

# updated cell annotation plot (Fig. 4E)
DimPlot(pbmc, group.by = "New_Annotations", label = F)

# make a heatmap of ADT markers
protein_marks <- c("CITE-CD3-", "CITE-CD45", "CITE-CD57-",
                   "CITE-CD27-", "CITE-CD28-", "CITE-CD161-",
                   "CITE-CD196", "CITE-CD127-", "CITE-CD183", "CITE-CD185",
                   "CITE-CD195", "CITE-CD194", "CITE-CD197", "CITE-CD25-", "CITE-CD69-",
                   "CITE-CD278")
matches <- unique (grep(paste(protein_marks,collapse="|"), 
                        rownames(pbmc), value=TRUE))

dp_genes <- DotPlot(pbmc, features = matches, group.by = "New_Annotations") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dp_genes
dp_genes_mat <- dp_genes$data
dp_genes_mat[1:4,]

# remove barcodes from cite-seq labels
dp_genes_mat$short <- substr(x = as.character(dp_genes_mat$features.plot), start = 1, stop = nchar(as.character(dp_genes_mat$features.plot))-16)
dp_genes_mat$short <- gsub(x = dp_genes_mat$short, pattern = "CITE-", replacement = "")

dp_genes_mat_sel <- dp_genes_mat[, c("short", "id", "avg.exp.scaled")]
dp_genes_mat_sel[1:4,]
rownames(dp_genes_mat_sel) <- 1:nrow(dp_genes_mat_sel)
exp_mat <- spread(dp_genes_mat_sel, key = id, avg.exp.scaled)

rownames(exp_mat) <- exp_mat$short
exp_mat[,1] <- NULL
dim(exp_mat)
exp_mat[1:4,]

# make a heatmap
pheatmap(mat = t(exp_mat), color = viridis(80),
         show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T,
         clustering_method = "ward.D2")

