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

# infer trajectory
model <- infer_trajectory(dataset, method = methods_selected[1], seed = 1, give_priors = c("start_id", "start_n", "end_n"))

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
colors_ann <- list(Treatment = c(Baseline = "#2657de", LPS = "#de2d26"),
                   Pseudotime = viridis(80))
anno_col <- anno_col[, 1:2]
anno_col <- anno_col[names(pseudo_ord), ]
breaksList = seq(-2, 2, by = 0.1)
length(breaksList)

# heatmap of genes changing with pseudotime (Fig. S2F)
pheatmap(mat = norm_sel[unique(asd$feature_id[1:100]), rownames(anno_col)], annotation_colors = colors_ann,
         show_rownames = T, show_colnames = F, scale = "row", cluster_cols = F,
         annotation_col = anno_col, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                              "RdBu")))(length(breaksList)),
         breaks = breaksList)


# histogram of cells along pseudotime (Fig. S2E)
anno_col$Treatment <- factor(anno_col$Treatment, levels = c("Baseline", "LPS"))
gp <- ggplot(anno_col, aes(x = Pseudotime, fill = Treatment)) +
  geom_histogram(bins = 10) +
  facet_wrap(~Treatment, ncol = 1) + 
  scale_fill_manual(values = c("#2657de" , "#de2d26")) +
  labs(y = "Number of Cells")
plot(gp)


# plot pseudotime of monocytes by donor (violin plots)
anno_col[1:4,]
pbmc_meta <- dataset$cell_info
table(rownames(anno_col) == rownames(pbmc_meta))
pbmc_meta_ord <- pbmc_meta[rownames(anno_col), ]
table(rownames(anno_col) == rownames(pbmc_meta_ord))
# add donor info
sp1 <- strsplit(x = pbmc_meta_ord$BEST.GUESS, split = ",")
sp2 <- sapply(sp1, `[[`, 1)
anno_col$Donor <- sp2

# plot pseudotime values of cells across each donor
gp <- ggplot(anno_col, aes(x = Donor, y = Pseudotime, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#2657de" , "#de2d26")) +
  labs(y = "Pseudotime", x = "Donor") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(gp)

