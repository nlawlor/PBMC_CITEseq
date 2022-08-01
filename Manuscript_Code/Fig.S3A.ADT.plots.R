# Re plot dim plots
library(Seurat)
library(ggplot2)
library(Matrix)
library(RColorBrewer)

# load the processed dataset
ser <- readRDS("../../CZI_Fourth_Run_no_doublets_pairwise.LPS.seurat.output.Rds")

# plot the cell surface proteins
id1 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD14-"))
id2 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD19-"))
id3 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD56-"))
id4 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD3-"))

# plot of cell surface markers
FeaturePlot(ser, features = rownames(ser)[c(id1, id2, id3, id4)], min.cutoff = "q05", max.cutoff = "q95")

id1 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD4-"))
id2 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD8a-"))
id3 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD45RA-"))
id4 <- which(grepl(x = rownames(ser), pattern = "^CITE-CD45RO-"))
FeaturePlot(ser, features = rownames(ser)[c(id1, id2, id3, id4)], min.cutoff = "q05", max.cutoff = "q95")
