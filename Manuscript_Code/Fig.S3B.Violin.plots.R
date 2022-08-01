# Monocyte analyze baseline and LPS cells, no doublets
rm(list = ls())
library(Matrix)
library(pheatmap)
library(Biobase)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(Vennerable)

ser <- readRDS("../../Monocyte.doublets.removed.seurat.output.Rds")
dim(ser)
meta_ser <- ser@meta.data
meta_ser <- meta_ser[meta_ser$Treatment %in% c("Control", "Mono"), ]
ser_sub <- subset(ser, cells = rownames(meta_ser))

counts_sel <- GetAssayData(ser_sub, slot = "counts")
meta_new <- ser_sub@meta.data
meta_use <- meta_new[, c("Treatment", "Cell_Type", "BEST.GUESS", "NUM.SNPS", "NUM.READS", "DROPLET.TYPE", "orig.ident")]
id_adt <- which(grepl(x = rownames(counts_sel), pattern = "^CITE-"))

pbmc <- CreateSeuratObject(counts = counts_sel, meta.data = meta_use, assay = "RNA", min.features = 0)
print("Creating seurat object")
pbmc <- NormalizeData(object = pbmc, normalization.method = "CLR")
# variable features
pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'vst', nfeatures = 500)
print(length(x = VariableFeatures(object = pbmc)))
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc))
print("Normalizing and scaling")
pbmc <- RunPCA(object = pbmc, features = union(rownames(counts_sel)[id_adt], VariableFeatures(object = pbmc)), verbose = FALSE)
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
# For Figure 2C
print(DimPlot(object = pbmc, reduction = 'umap', group.by = "Treatment", cols = c("#2657de", "#de2d26")))
print(DimPlot(object = pbmc, reduction = 'umap', group.by = "Cell_Type"))


# additional feature plots (for Fig S2C)
VlnPlot(pbmc, features = c("IL8", "IL1B"), group.by = "Treatment", cols = c("#2657de", "#de2d26"))

