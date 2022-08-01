# Re plot dim plots
library(Seurat)
library(ggplot2)
library(Matrix)
library(RColorBrewer)

cols <- c("#fdbf6f", "#ff7f00", "#f7fcb9", "#addd8e", "#41ab5d", "#006837",
          "lightskyblue1", "lightskyblue3", "#2171b5", "#08306b", 
          "rosybrown1", "red2", "red4", "#c2a5cf", "#762a83")
celltype <- c("B", "Actv_B", "CD4T_Naive", "Actv_CD4T_Naive", "CD4T_Mem", "Actv_CD4T_Mem",
              "CD8T_Naive", "Actv_CD8T_Naive", "CD8T_Mem", "Actv_CD8T_Mem",
              "CD16_Mono", "CD14_Mono", "Actv_CD14_Mono", "NK", "Actv_NK")

# do for mono
ser <- readRDS("../../CZI_Fourth_Run_no_doublets_pairwise.LPS.seurat.output.Rds")
meta <- ser@meta.data

# identify which colors
id_c <- which(celltype %in% meta$Annotation)

# plot dimplot
ser_dim <- as.data.frame(ser@reductions$umap@cell.embeddings)
table(rownames(ser_dim) == colnames(ser))
ser_dim$Annotation <- meta$Annotation
ser_dim <- within(ser_dim,
                  Annotation <- factor(Annotation,
                                     levels=celltype[id_c]))

# plot of cell annotations
p1 <- ggplot(ser_dim, aes(x = UMAP_1, y = UMAP_2, color = Annotation)) +
  geom_point(size = 1, alpha = 0.5) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        axis.title=element_text(size=14,face="bold")) +
  theme(legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values = cols[id_c], labels = celltype[id_c]) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "UMAP 1", y = "UMAP 2")
plot(p1)


# plot by treatments
ser_dim$Treatment <- meta$Treatment
p1 <- ggplot(ser_dim, aes(x = UMAP_1, y = UMAP_2, color = Treatment)) +
  geom_point(size = 1, alpha = 0.5) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        axis.title=element_text(size=14,face="bold")) +
  theme(legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values = c("#2657de", "#de2d26")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "UMAP 1", y = "UMAP 2")
plot(p1)