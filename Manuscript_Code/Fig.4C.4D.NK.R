rm(list = ls())
set.seed(100)
library(Rtsne)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Matrix)

# read in old seurat file
ser <- readRDS("Nk.cells.no.doublets.no.IgM.seurat.Rds")
dim(ser)
meta <- ser@meta.data
meta_sel <- meta[meta$Treatment %in% c("Control", "Tcell"), ]
dim(meta_sel)
ser_sub <- subset(ser, cells = rownames(meta_sel))
meta_sel <- meta_sel[, c("orig.ident", "Treatment", "Cell_Type",
                         "BEST.GUESS", "NUM.SNPS", "NUM.READS",
                         "DROPLET.TYPE")]

# get counts
comb.counts <- GetAssayData(ser_sub, slot = "counts")
dim(comb.counts)

idx <- which(grepl(x = rownames(comb.counts), pattern = "^CITE-"))
adts <- rownames(comb.counts)[idx]

pdf(file = "NK.baseline.antiCD3.plots.pdf", onefile = T)
# create object and split split the data by treatment type
pbmc <- CreateSeuratObject(counts = comb.counts, meta.data = meta_sel, assay = "RNA", min.features = 0)
print("Creating seurat object")
pbmc <- NormalizeData(object = pbmc, normalization.method = "CLR")
# variable features
pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'vst', nfeatures = 500)
print(length(x = VariableFeatures(object = pbmc)))
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc))
print("Normalizing and scaling")
pbmc <- RunPCA(object = pbmc, features = union(adts, VariableFeatures(object = pbmc)), verbose = FALSE)
#print(DimPlot(object = pbmc))
print("Running PCA")

pbmc <- RunUMAP(object = pbmc, dims = 1:20)
pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
pbmc <- FindClusters(object = pbmc, resolution = 1.2)

# make a histogram of number of features
print(hist(as.numeric(pbmc$nFeature_RNA), xlab = "Number of Genes Detected", main = ""))
summ <- summary(pbmc$nFeature_RNA)
legend("topright", legend = paste(names(summ), summ, sep = " "))

# categorize detect gene info
det_gen <- pbmc$nFeature_RNA
Detect_Gene <- cut(det_gen, breaks = 20)
pbmc <- AddMetaData(pbmc, metadata = Detect_Gene, col.name = "Detected_Genes")

print(DimPlot(object = pbmc, reduction = 'umap', label = TRUE))
print(DimPlot(object = pbmc, reduction = 'umap', group.by = "Treatment"))
print(DimPlot(object = pbmc, reduction = 'umap', group.by = "Treatment", 
              cols = c("#377eb8", "#984ea3")))
print(DimPlot(object = pbmc, reduction = 'umap', group.by = "Cell_Type", label = TRUE))
#print(DimPlot(object = pbmc, reduction = 'umap', group.by = "Annotation", label = TRUE))

# print marker genes on tsne plot
marker_genes <- c("NKG7", "GZMH", "CCL5", "PRF1", "GNLY", "KLRB1", "MZB1", 
                  "IGLL5", "CD8A", "CD8B", "CD4", "CD27", "CD3", "CD3D", "CD3E", "IL7R", "CD56", "FOXP3", "CD25",
                  "CD14", "CD16", "LYZ", "CST3", "FCER1A",
                  "S100A8", "S100A9", "CD19", "MS4A1", "CD79A", "IGJ", "FCGR3A", "MS4A7",
                  "FCER1G", "LILRA4", "IRF7",
                  "PPBP", "PF4", "HBA2", "HBA1", "HBB", "CD34",
                  "CD1C", "CD141")

# print features
all_feats <- c(marker_genes, adts)
id_feat <- which(all_feats %in% rownames(pbmc))
all_feats_sel <- all_feats[id_feat]
for (f in seq(from = 1, to = length(all_feats_sel), by = 4)) {
  print(FeaturePlot(pbmc, features = all_feats_sel[(f):(f+3)], min.cutoff = "q05", max.cutoff = "q95", pt.size = 0.15))
  #print(f:(f+3))
}

# find rna markers and make a heatmap
# rna.markers <- FindAllMarkers(pbmc, only.pos = F, 
#                               print.bar = F)
# colnames(rna.markers)[7] <- "Marker"
# write.csv(x = rna.markers, file = paste(name, treats[t], "markers.csv", sep = "."))

# cell proportions by cluster
# tab <- table(pbmc$RNA_snn_res.1.2, pbmc$Treatment)
# print(tab)
# write.csv(x = tab, file = paste(name, treats[t], "cell.proportions.csv", sep = "."))

saveRDS(pbmc, file = paste("NK.baseline.antiCD3.seurat.output.Rds", sep = "."))

dev.off()

pdf(file = "NK.baseline.antiCD3.dimplot.pdf", onefile = T)
DimPlot(object = pbmc, reduction = 'umap', group.by = "Treatment", 
              cols = c("#377eb8", "#e41a1c"), pt.size = 2)
DimPlot(object = pbmc, reduction = 'umap', pt.size = 2)
dev.off()

tab_clus <- as.data.frame(table(pbmc@meta.data$seurat_clusters, pbmc@meta.data$Treatment))
tab_clus[1:4,]

pdf(file = "NK.baseline.antiCD3.clusters.pdf")
ggplot(tab_clus, aes(x=Var1, y=Freq, fill=Var2)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
  labs(x = "Cluster", y = "% of cells per treatment", fill = "Treatment")
dev.off()


# heatmap of genes in order, no row clustering
genes <- rev(c("IFI6", "IFI27", "LY6E", "ISG15", "IFI44L", "STAT1", "MX1",
           "IFITM1", "IFIT3", "CD69", "CD83", "XCL1", "XCL2", 
           "CCL3", "CCL4", "KLF2", "KLF4", "FCRL6", "FGFBP2",
           "CD52", "CCL5"))
dp <- DotPlot(pbmc, features = genes)
# change from long to wide and make pheatmap
long_df <- dp$data
library(tidyr)
data_wide <- spread(long_df[, c("avg.exp", "id", "features.plot")], id, avg.exp)
dim(data_wide)
data_wide[1:4,]
rownames(data_wide) <- data_wide[,1]
data_wide[,1] <- NULL
pdf(file = "NK.baseline.antiCD3.heatmap.pdf", onefile = T)
pheatmap(data_wide, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                             "RdBu")))(100),
         cluster_rows = F, cluster_cols = T,
         clustering_method = "ward.D2",
         scale = "row")
dev.off()


