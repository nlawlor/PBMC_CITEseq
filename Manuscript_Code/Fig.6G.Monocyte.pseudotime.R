rm(list = ls())
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
pbmc <- readRDS(file = "../../Monocyte.doublets.removed.seurat.output.Rds")

# get select cells
pbmc_meta <- pbmc@meta.data
dim(pbmc_meta)
pbmc_meta <- pbmc_meta[pbmc_meta$Treatment %in% c("Control", "Mono"), ]
dim(pbmc_meta)
# add editied treatment info
pbmc_meta$Stimulation <- pbmc_meta$Treatment
pbmc_meta$Stimulation <- gsub(x = pbmc_meta$Stimulation, pattern = "Mono", replacement = "LPS")
pbmc_meta$Stimulation <- gsub(x = pbmc_meta$Stimulation, pattern = "Control", replacement = "Baseline")
pbmc_sub <- subset(pbmc, cells = rownames(pbmc_meta))

counts <- GetAssayData(pbmc_sub, slot = "counts")
norm_data <- GetAssayData(pbmc_sub, slot = "data")


mono_all <- read.csv("../../Monocyte.doublets.removed.LPS.Induced.unique.markers.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)

# get genes induced with treatment
mono_genes <- rownames(mono_all[mono_all$avg_logFC > 0.25 & mono_all$p_val_adj < 0.05, ])

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

guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected
methods_selected

model <- infer_trajectory(dataset, method = methods_selected[1], seed = 1, give_priors = c("start_id", "start_n", "end_n"))


# for Fig 2G (top) panel
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_density = "grouping",
  grouping = dataset$cell_info$Stimulation,
  label_milestones = F,
  alpha_cells = 0.5
) + scale_color_manual(values = c("#2657de", "#de2d26")) + scale_fill_manual(values = c("#2657de", "#de2d26"))

# calculate pseudotime
pseudo <- calculate_pseudotime(model)
plot_dimred(model, "pseudotime", pseudotime = calculate_pseudotime(model)) + ggtitle("Pseudotime")
# change inf to zero
id_inf <- which(is.infinite(pseudo))
length(id_inf)
pseudo[id_inf] <- 0
pseudo_ord <- sort(pseudo, decreasing = F)

# gene expression plots Fig 2G (bottom)
patchwork::wrap_plots(
  plot_dimred(model, feature_oi = "S100A9", expression_source = dataset) + ggtitle("S100A9"),
  plot_dimred(model, feature_oi = "S100A8", expression_source = dataset) + ggtitle("S100A8"),
  plot_dimred(model, feature_oi = "MT2A", expression_source = dataset) + ggtitle("MT1F")
)

patchwork::wrap_plots(
  plot_dimred(model, feature_oi = "IL1B", expression_source = dataset) + ggtitle("IL1B"),
  plot_dimred(model, feature_oi = "IL8", expression_source = dataset) + ggtitle("IL8"),
  plot_dimred(model, feature_oi = "IL6", expression_source = dataset) + ggtitle("IL6")
)
