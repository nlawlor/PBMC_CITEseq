# UpsetR plot of all response genes
library(ggplot2)
library(VennDiagram)
library(Vennerable)
library(UpSetR)

# load in gene lists 
nk_cd3 <- read.csv("NK.doublets.removed.AntiCD3.vs.Unactivated.marker.genes.csv",
                   header = T, check.names = F, stringsAsFactors = F, row.names = 1)
nk_cd3_up <- nk_cd3[nk_cd3$avg_logFC > 0.25 & nk_cd3$p_val_adj < 0.05, ]

cd4t_naive <- read.csv("CD4T.doublets.removed.Activated.Naive.vs.Unactivated.marker.genes.csv",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
cd4t_naive_up <- cd4t_naive[cd4t_naive$avg_logFC > 0.25 & cd4t_naive$p_val_adj < 0.05, ]

cd4t_mem <- read.csv("CD4T.doublets.removed.Activated.Memory.vs.Unactivated.marker.genes.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
cd4t_mem_up <- cd4t_mem[cd4t_mem$avg_logFC > 0.25 & cd4t_mem$p_val_adj < 0.05, ]

cd8t_naive <- read.csv("CD8T.doublets.removed.Activated.Naive.vs.Unactivated.marker.genes.csv",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
cd8t_naive_up <- cd8t_naive[cd8t_naive$avg_logFC > 0.25 & cd8t_naive$p_val_adj < 0.05, ]

cd8t_mem <- read.csv("CD8T.doublets.removed.Activated.Memory.vs.Unactivated.marker.genes.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
cd8t_mem_up <- cd8t_mem[cd8t_mem$avg_logFC > 0.25 & cd8t_mem$p_val_adj < 0.05, ]

# create a list of the induced genes for NK and T cells
induced_genes <- list(rownames(nk_cd3_up),
                   rownames(cd4t_naive_up),
                   rownames(cd4t_mem_up),
                   rownames(cd8t_naive_up),
                   rownames(cd8t_mem_up))
names(induced_genes) <- c("NK", 
                 "CD4T_Naive", "CD4T_Mem", "CD8T_Naive", "CD8T_Mem")

upset(fromList(induced_genes), order.by = "freq",
      nsets = length(induced_genes), nintersects = 20,
      mainbar.y.label = "# of Gene Intersections",
      sets.x.label = "# of Induced Genes", mb.ratio = c(0.55, 0.45),
      text.scale = c(1.5, 2, 1.2, 1.2, 1.5, 1.5))

