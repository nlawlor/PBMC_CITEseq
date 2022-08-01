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

# load seurat object
ser <- readRDS("../../CZI_Fourth_Run_no_doublets_pairwise.antiCD3CD28.seurat.output.Rds")
meta <- ser@meta.data

# identify which colors
id_c <- which(celltype %in% meta$Annotation)

ser_dim <- as.data.frame(ser@reductions$umap@cell.embeddings)

# dimension reduction plot by cell annotation, for Fig. 4A
p1 <- ggplot(ser_dim, aes(x = UMAP_1, y = UMAP_2, color = Annotation)) +
  geom_point(size = 0.15, alpha = 0.5) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        axis.title=element_text(size=14,face="bold")) +
  theme(legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values = cols[id_c], labels = celltype[id_c]) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "UMAP 1", y = "UMAP 2")
plot(p1)



# plot dimensions and color by treatment
ser_dim$Treatment <- meta$Treatment
p1 <- ggplot(ser_dim, aes(x = UMAP_1, y = UMAP_2, color = Treatment)) +
  geom_point(size = 0.15, alpha = 0.5) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        axis.title=element_text(size=14,face="bold")) +
  theme(legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values = c("#de2d26", "#2657de")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "UMAP 1", y = "UMAP 2")
plot(p1)
dev.off()


# plot cite-seq markers for Fig 4B
id1 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD4-"))
id2 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD8a-"))
id3 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD3-"))

print(FeaturePlot(ser, features = rownames(ser)[c(id1, id2, id3)],
                  cols = c("lightgrey", "blue"),
                  min.cutoff = "q05", max.cutoff = "q95"))


id1 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD45RA-"))
id2 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD45RO-"))
id3 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD25-"))
id4 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD69-"))

print(FeaturePlot(ser, features = rownames(ser)[c(id1, id2, id3, id4)],
                  cols = c("lightgrey", "blue"),
                  min.cutoff = "q05", max.cutoff = "q95"))
