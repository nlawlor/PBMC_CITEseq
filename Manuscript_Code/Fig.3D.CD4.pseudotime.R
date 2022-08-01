library(Matrix)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dyno)
library(tidyverse)
library(dplyr)
library(viridis)

# for cd8 t cells first
pbmc <- readRDS("../../CD4T.doublets.removed.no.IgM.seurat.Rds")

# get select cells
pbmc_meta <- pbmc@meta.data
dim(pbmc_meta)
counts <- GetAssayData(pbmc, slot = "counts")
norm_data <- GetAssayData(pbmc, slot = "data")

# get only baseline and stim
pbmc_meta <- pbmc_meta[pbmc_meta$Treatment %in% c("Baseline", "Anti_CD3_CD28"), ]

# load in differential genes
mem_cd3 <- read.csv("../..//CD4T.doublets.removed.Activated.Memory.vs.Unactivated.marker.genes.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
mem_cd3_up <- mem_cd3[mem_cd3$avg_logFC > 0.25 & mem_cd3$p_val_adj < 0.05,]

naive_cd3 <- read.csv("../../CD4T.doublets.removed.Activated.Naive.vs.Unactivated.marker.genes.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
naive_cd3_up <- naive_cd3[naive_cd3$avg_logFC > 0.25 & naive_cd3$p_val_adj < 0.05,]

# take union of induced genes
up_genes <- union(rownames(mem_cd3_up), rownames(naive_cd3_up))

# subset matrices by these genes
counts_sel <- counts[up_genes, rownames(pbmc_meta)]
norm_sel <- norm_data[up_genes, rownames(pbmc_meta)]
dim(counts_sel)
pbmc_meta$cell_id <- rownames(pbmc_meta)

# create dataset
dataset <- wrap_expression(
  counts = t(counts_sel),
  expression = t(norm_sel),
  cell_info = pbmc_meta
)

# add prior info (control cells)
pbmc_cont_naive <- dataset$cell_info[dataset$cell_info$Treatment == "Baseline" & dataset$cell_info$New_Annotations == "CD4T_Naive",]
pbmc_cont_mem <- dataset$cell_info[dataset$cell_info$Treatment == "Baseline" & dataset$cell_info$New_Annotations == "CD4T_Mem",]
pbmc_treat_naive <- dataset$cell_info[dataset$cell_info$New_Annotations == "CD4T_Naive_Actv", ]
pbmc_treat_mem <- dataset$cell_info[dataset$cell_info$New_Annotations == "CD4T_Mem_Actv", ]

dataset <- add_prior_information(dataset = dataset, start_id = c(sample(rownames(pbmc_cont_naive), 1), sample(rownames(pbmc_cont_mem), 1)),
                                 start_n = 2, end_n = 2, end_id = c(sample(rownames(pbmc_treat_naive), 1), sample(rownames(pbmc_treat_mem), 1)))
guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected
methods_selected
model <- infer_trajectory(dataset, method = methods_selected[1], seed = 1,
                          give_priors = c("start_id", "start_n", "end_n", "end_id"))

# plot trajectory
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_density = "none",
  grouping = dataset$cell_info$New_Annotations,
  label_milestones = F,
  alpha_cells = 0.5
) 

plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_density = "none",
  grouping = dataset$cell_info$Treatment,
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

anno_col <- data.frame(Treatment = dataset$cell_info$Treatment,
                       Annotation = dataset$cell_info$New_Annotations,
                       Pseudotime = pseudo,
                       Cell_ID = dataset$cell_info$cell_id)
rownames(anno_col) <- anno_col$Cell_ID
anno_col <- anno_col[, 1:3]
anno_col <- anno_col[names(pseudo_ord), ]
anno_col <- anno_col[, c("Treatment", "Pseudotime", "Annotation")]


# plot pseudotime of monocytes by donor (violin plots)
anno_col[1:4,]
pbmc_meta_ord <- pbmc_meta[rownames(anno_col), ]
table(rownames(anno_col) == rownames(pbmc_meta_ord))
# add donor info
sp1 <- strsplit(x = pbmc_meta_ord$BEST.GUESS, split = ",")
sp2 <- sapply(sp1, `[[`, 1)
anno_col$Donor <- sp2
gp <- ggplot(anno_col, aes(x = Annotation, y = Pseudotime, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#984ea3", "grey" )) +
  labs(y = "Pseudotime", x = "Annotation") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~Annotation)
plot(gp)

