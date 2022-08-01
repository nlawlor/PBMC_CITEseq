rm(list = ls())
library(Seurat)
library(UpSetR)
library(Vennerable)
library(ggplot2)
library(dyno)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(egg)

# load model and dataset
model <- readRDS("../../Monocyte.LPS.and.control.induced.genes.pseudotime.slingshot.model.Rds")
dataset <- readRDS("../../Monocyte.LPS.and.control.induced.genes.pseudotime.dataset.Rds")

# calculate pseudotime
pseudo <- calculate_pseudotime(model)
# change inf to zero
id_inf <- which(is.infinite(pseudo))
length(id_inf)
pseudo[id_inf] <- 0
pseudo_ord <- sort(pseudo, decreasing = F)

# data frame of cell information
anno_col <- data.frame(Treatment = dataset$cell_info$Stimulation,
                       Pseudotime = pseudo,
                       Cell_ID = dataset$cell_info$cell_id)

# do kmeans clustering to divide cells into 4 groups
set.seed(1)
km <- kmeans(model$dimred[, 1:2], centers = 4, iter.max = 1000, nstart = 1000,
             algorithm = "Hartigan-Wong")
mydata <- data.frame(model$dimred, km$cluster)
table(rownames(mydata) == rownames(anno_col))
km_ord <- mydata$km.cluster
id1 <- which(km_ord == 1)
km_ord[id1] <- "LPS_3"
id4 <- which(km_ord == 4)
km_ord[id4] <- "Baseline"
id2 <- which(km_ord == 3)
km_ord[id2] <- "LPS_1"
id3 <- which(km_ord == 2)
km_ord[id3] <- "LPS_2"
anno_col$Kmeans <- km_ord
table(anno_col$Kmeans)


grp1 <- which(anno_col$Treatment == "LPS" & anno_col$Pseudotime > 5 & anno_col$Pseudotime < 13)
grp2 <- which(anno_col$Treatment == "LPS" & anno_col$Pseudotime > 13 & anno_col$Pseudotime < max(anno_col$Pseudotime))
grp3 <- which(anno_col$Treatment == "LPS" & anno_col$Pseudotime == max(anno_col$Pseudotime))
anno_col$Cluster <- "Baseline"
anno_col$Cluster[grp1] <- "LPS_1"
anno_col$Cluster[grp2] <- "LPS_2"
anno_col$Cluster[grp3] <- "LPS_3"
table(anno_col$Cluster)

# Fig 3A
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_density = "grouping",
  grouping = anno_col$Cluster,
  label_milestones = F,
  alpha_cells = 0.5
) + ggtitle("Monocyte States") + 
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + 
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))

# read in seurat object
ser <- readRDS("../../LPS.vs.Baseline.monocytes.no.doublets.seurat.Rds")
dim(ser)
table(colnames(ser) == anno_col$Cell_ID)
ser <- AddMetaData(ser, col.name = "Pseudotime_Cluster", metadata = anno_col$Cluster)
ser <- AddMetaData(ser, col.name = "Pseudotime", metadata = anno_col$Pseudotime)
ser <- AddMetaData(ser, col.name = "Stimulation", metadata = anno_col$Treatment)

# read in inflammation scores
inflamm_scores <- readRDS("../../LPS.vs.Baseline.monocytes.no.doublets.seurat.inflamm.scores.Rds")
dim(inflamm_scores)
table(colnames(ser) == rownames(inflamm_scores))
ser <- AddMetaData(ser, col.name = "Inflamm_Score", metadata = inflamm_scores$Inflamm_Score)

# plot of inflammation scores; Fig 3C
RidgePlot(ser, features = "Inflamm_Score", group.by = "Pseudotime_Cluster")

# barplot of the number of cells in each state, also proportions per donor
table(ser$Pseudotime_Cluster)
Donor <- as.character(ser$BEST.GUESS)
d1 <- strsplit(x = Donor, split = "\\,")
Donor_Id <- sapply(d1, `[[`, 1)
ser <- AddMetaData(ser, col.name = "Donor_Id", metadata = Donor_Id)
tab_cell_num <- as.data.frame(table(ser$Pseudotime_Cluster))

# Fig S3A, plot of cell numbers in each state
ggplot(tab_cell_num, aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  geom_text(aes(y = Freq, label = Freq)) + 
  labs(y = "# of Monocytes in each State", x = "")

# Plot pseudotime values of cells in each state; Fig S3B
VlnPlot(ser, features = "Pseudotime", group.by = "Pseudotime_Cluster") + coord_flip()

# Fig 3B (plot of cell proportions in each state)
tab_don <- as.data.frame(table(ser$Pseudotime_Cluster, ser$Donor_Id))
ggplot(tab_don, aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "# of Monocytes in each State per Donor", x = "Donor of Origin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# use the semi-supervised clustering results identify marker genes
Idents(ser) <- ser$Pseudotime_Cluster
marks_all <- FindAllMarkers(ser, only.pos = T, logfc.threshold = 0.25)
dim(marks_all)
table(marks_all$cluster)

# make upset plot of genes
id_cite <- which(grepl(x = marks_all$gene, pattern = "^CITE-"))
mark1 <- marks_all[-id_cite, ]

listInput <- list()
for (l in 1:length(unique(mark1$cluster))) {
  # get genes 
  sel <- mark1[mark1$cluster == as.character(unique(mark1$cluster)[l]),]
  listInput[[l]] <- sel$gene
  names(listInput)[l] <- as.character(unique(mark1$cluster)[l])
}

# Fig. S3C, upset diagram of genes induced in each state
upset(fromList(listInput), order.by = "freq",
      nsets = length(listInput), nintersects = 20,
      mainbar.y.label = "# of Gene Intersections",
      sets.x.label = "# of Induced Genes", mb.ratio = c(0.55, 0.45),
      text.scale = c(1.5, 2, 1.2, 1.2, 1.5, 1))

# get genes unique to each state(s)
venn.snp <- Venn(SetNames = names(listInput),
                 Sets =listInput)
sapply(venn.snp@IntersectionSets, length)
LPS_3_only <- venn.snp@IntersectionSets$`1000`
LPS_2_only <- venn.snp@IntersectionSets$`0100`
baseline_only <- venn.snp@IntersectionSets$`0001`
LPS_2_3 <- venn.snp@IntersectionSets$`1100`


# re-annotate the lps2/3 common genes (if greater than logfc 0.5 difference, put in unique group)
lps3_new <- NULL
lps_2_new <- NULL
for (g in 1:length(LPS_2_3)) {
  md3 <- mark1[mark1$gene == LPS_2_3[g] & mark1$cluster == "LPS_3",]
  md2 <- mark1[mark1$gene == LPS_2_3[g] & mark1$cluster == "LPS_2",]
  if (md3$avg_logFC - md2$avg_logFC >= 0.5) {
    lps3_new <- c(lps3_new, LPS_2_3[g])
  } else if (md2$avg_logFC - md3$avg_logFC >= 0.5) {
    lps_2_new <- c(lps_2_new, LPS_2_3[g])
  } else {}
}

id_rem <- which(LPS_2_3 %in% c(lps3_new, lps_2_new))
LPS_2_3 <- LPS_2_3[-id_rem]
LPS_3_only <- c(LPS_3_only, lps3_new)
LPS_2_only <- c(LPS_2_only, lps_2_new)


meta_df <- ser@meta.data
# get selected metadata for plot
meta_sel <- meta_df[, c("Stimulation", "Pseudotime_Cluster")]

# get norm data
norm_data <- GetAssayData(ser, slot = "data")
dim(norm_data)

# heatmap of genes by pseudotime
c1 <- which(meta_sel$Pseudotime_Cluster == "Baseline")
c2 <- which(meta_sel$Pseudotime_Cluster == "LPS_1")
c3 <- which(meta_sel$Pseudotime_Cluster == "LPS_2")
c4 <- which(meta_sel$Pseudotime_Cluster == "LPS_3")
  
colors_ann <- list(Stimulation = c(Baseline = "#2657de", LPS = "#de2d26"))


# make plots of the induced genes for each state
# state LPS3 genes
dp3 <- DotPlot(ser, features = sort(LPS_3_only), group.by = "Pseudotime_Cluster",
               cols = c("lightgrey", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 8),
        axis.title.y = element_blank()) + coord_flip() + NoLegend() +
  ggtitle("LPS_3 Induced")

# make a plot with legend to extract
dp3_with_leg <- DotPlot(ser, features = sort(LPS_3_only), group.by = "Pseudotime_Cluster",
                        cols = c("lightgrey", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.position = "bottom") + coord_flip() +
  ggtitle("LPS_3 Induced")
# extract legend to reuse later
leg <- cowplot::get_legend(dp3_with_leg)


# LPS2 state genes
dp2 <- DotPlot(ser, features = sort(LPS_2_only), group.by = "Pseudotime_Cluster",
               cols = c("lightgrey", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 8),
        axis.title.y = element_blank()) + coord_flip() +
  ggtitle("LPS_2 Induced") + NoLegend()

# LPS2 and LPS3 genes
dp2_3 <- DotPlot(ser, features = sort(LPS_2_3), group.by = "Pseudotime_Cluster",
                 cols = c("lightgrey", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 8),
        axis.title.y = element_blank()) + coord_flip() +
  ggtitle("LPS_2 and LPS_3 Induced") + NoLegend()


# make heatmaps of scaled expression
dp3_sel <- dp3$data
dp3_sel <- dp3_sel[, c("id", "features.plot",  "avg.exp.scaled")]
dp3_df <- tidyr::spread(dp3_sel, id, avg.exp.scaled)
rownames(dp3_df) <- dp3_df$features.plot
dp3_df[,1] <- NULL
ph3 <- pheatmap(dp3_df, cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                           "RdBu")))(100),
         show_rownames = T, show_colnames = T)

dp2_sel <- dp2$data
dp2_sel <- dp2_sel[, c("id", "features.plot",  "avg.exp.scaled")]
dp2_df <- tidyr::spread(dp2_sel, id, avg.exp.scaled)
rownames(dp2_df) <- dp2_df$features.plot
dp2_df[,1] <- NULL
ph2 <- pheatmap(dp2_df, cluster_rows = T, cluster_cols = F,
                color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                          "RdBu")))(100),
                show_rownames = T, show_colnames = T)

dp2and3_sel <- dp2_3$data
dp2and3_sel <- dp2and3_sel[, c("id", "features.plot",  "avg.exp.scaled")]
dp2and3_df <- tidyr::spread(dp2and3_sel, id, avg.exp.scaled)
rownames(dp2and3_df) <- dp2and3_df$features.plot
dp2and3_df[,1] <- NULL
ph2and3 <- pheatmap(dp2and3_df, cluster_rows = T, cluster_cols = F,
                color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                          "RdBu")))(100),
                show_rownames = T, show_colnames = T)


# combined  heatmap of induced genes; Fig S3D
g <- grid.arrange(arrangeGrob(grobs= list(ph2[[4]], ph2and3[[4]],
                                          ph3[[4]]),ncol=3))

