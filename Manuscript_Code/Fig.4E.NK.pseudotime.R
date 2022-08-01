library(Matrix)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dyno)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(viridis)

# load in object
pbmc <- readRDS(file = "../../Nk.cells.no.doublets.no.IgM.seurat.Rds")

# get select cells
pbmc_meta <- pbmc@meta.data
pbmc_meta <- pbmc_meta[pbmc_meta$Treatment %in% c("Control", "Tcell"), ]
dim(pbmc_meta)
# add editied treatment info
pbmc_meta$Stimulation <- pbmc_meta$Treatment
pbmc_meta$Stimulation <- gsub(x = pbmc_meta$Stimulation, pattern = "Tcell", replacement = "Anti_CD3_CD28")
pbmc_meta$Stimulation <- gsub(x = pbmc_meta$Stimulation, pattern = "Control", replacement = "Baseline")
pbmc_sub <- subset(pbmc, cells = rownames(pbmc_meta))

counts <- GetAssayData(pbmc_sub, slot = "counts")
norm_data <- GetAssayData(pbmc_sub, slot = "data")

mono_all <- read.csv("../../NK.doublets.removed.AntiCD3.vs.Unactivated.marker.genes.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)

# load in differential genes
mono_up <- mono_all[mono_all$avg_logFC > 0.25 & mono_all$p_val_adj < 0.05,]
mono_genes <- rownames(mono_up)

counts_sel <- counts[mono_genes, rownames(pbmc_meta)]
norm_sel <- norm_data[mono_genes, rownames(pbmc_meta)]
dim(counts_sel)
pbmc_meta$cell_id <- rownames(pbmc_meta)

# create dataset
dataset <- wrap_expression(
  counts = t(counts_sel),
  expression = t(norm_sel),
  cell_info = pbmc_meta
)


# add prior info (control cells)
pbmc_cont <- dataset$cell_info[dataset$cell_info$Treatment == "Control",]
dim(pbmc_cont)
dataset <- add_prior_information(dataset = dataset, start_id = rownames(pbmc_cont), start_n = 1,
                                 end_n = 1)

# choose pseudotime method
guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected
methods_selected

# infer trajectory
model <- infer_trajectory(dataset, method = methods_selected[1], seed = 1, give_priors = c("start_id", "start_n", "end_n"))

# plot the trajectory
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_density = "grouping",
  grouping = dataset$cell_info$Stimulation,
  label_milestones = F,
  alpha_cells = 0.5
) + scale_color_manual(values = c("#984ea3", "grey")) + scale_fill_manual(values = c("#984ea3", "grey"))

# calculate pseudotime
pseudo <- calculate_pseudotime(model)
plot_dimred(model, "pseudotime", pseudotime = calculate_pseudotime(model)) + ggtitle("Pseudotime")
# change inf to zero
id_inf <- which(is.infinite(pseudo))
length(id_inf)
pseudo[id_inf] <- 0
pseudo_ord <- sort(pseudo, decreasing = F)

# heatmap of genes by pseudotime
asd = dynfeature::calculate_milestone_feature_importance(trajectory = model, expression_source = dataset)
anno_col <- data.frame(Treatment = dataset$cell_info$Stimulation,
                       Pseudotime = pseudo,
                       Cell_ID = dataset$cell_info$cell_id)
rownames(anno_col) <- anno_col$Cell_ID
colors_ann <- list(Treatment = c(Baseline = "grey", Anti_CD3_CD28 = "#984ea3"),
                   Pseudotime = viridis(80))
anno_col <- anno_col[, 1:2]
anno_col <- anno_col[names(pseudo_ord), ]
breaksList = seq(-2, 2, by = 0.1)
length(breaksList)

pheatmap(mat = norm_sel[unique(asd$feature_id[1:40]), rownames(anno_col)], annotation_colors = colors_ann,
         show_rownames = T, show_colnames = F, scale = "row", cluster_cols = F,
         annotation_col = anno_col, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                              "RdBu")))(length(breaksList)),
         breaks = breaksList)

