# load in CITE-seq data
library(data.table)
library(Matrix)
library(Seurat)
library(mclust)
rm(list = ls())
name <- "CZI.PMBC.HTO.pooled"
all_hto <- readRDS("CZI.PMBC.HTO.pooled.HTO.matrix.Rds")

# CLR normalize the data for GMM modeling
sr_obj <- CreateSeuratObject(counts = all_hto[c("control-GTCAACTCTTTAGCG", "Mono-AGTAAGTTCAGCGTA", "Bcell-TGATGGCCTATTGGG", "Tcell-TTCCGCCTCTCTTTG"),])
# centered log-normalization
sr_obj <- NormalizeData(sr_obj, normalization.method = "CLR") # norm data in @data slot
#norm_counts <- sr_obj@data
norm_counts <- GetAssayData(sr_obj, slot = "data")
rownames(norm_counts) <- c("Control", "Mono", "Bcell", "Tcell")

# use gmm method to classify cell types
p.exp <- norm_counts
p.bin <- p.exp
p.bin[] <- 0
p.names <- rownames(p.exp)

# keep list of output paramters
params <- list()
for(i in 1:nrow(p.exp))
{
  x <- Mclust(p.exp[i,],G=2)
  p.bin[i,] <- x$classification
  params[[i]] <- x$parameters
  names(params)[i] <- rownames(norm_counts)[i]
}


c.sum <- apply(p.bin,2,sum)
sample.type <- character(length = ncol(p.bin))
# 4 markers, if sum is 4, then sample cannot be classified
sample.type[c.sum == 4] <- "None"
for(i in 1:length(c.sum)) {
  if(c.sum[i] == 5 & length(which(p.bin[,i] ==2) == 1)) {
    sample.type[i] <- p.names[which(p.bin[,i] ==2)]
  }
  
  if(c.sum[i] > 5 & length(which(p.bin[,i] ==2) == 2)) {
    ids.sel <- p.names[which(p.bin[,i] ==2)]
    nams <- ""
    for (num in 1:length(ids.sel)) {
      nams <- paste(ids.sel[num], nams, sep="+")  
    }
    nams <- gsub(x = nams, pattern = "\\++$", replacement = "")
    sample.type[i] <- nams
  }
}

# short sample assignment
d_ids <- which(grepl(x = sample.type, pattern = "\\+"))
Sample_Short <- sample.type
Sample_Short[d_ids] <- "Doublet"

# data frame of cell barcodes and assignments
s_anns <- data.frame(Cell_Barcode = colnames(p.bin), Sample = sample.type, Sample_Short = Sample_Short)
table(s_anns$Sample)
table(s_anns$Sample_Short)
# add run identifiers
runs <- strsplit(x = as.character(s_anns$Cell_Barcode), split = "_")
run_nam <- unlist(lapply(runs, `[[`, 1))
s_anns$Run <- run_nam
# proportions of each sample type
tab_s <- table(s_anns$Run, s_anns$Sample_Short)
rsums <- rowSums(tab_s)
per_s <- 100 * (tab_s / rsums)

write.csv(x = table(s_anns$Run, s_anns$Sample_Short), file = paste(name, "sample.classifications.GMM.summary.csv", sep = "."))
write.csv(x = per_s, file = paste(name, "sample.classifications.GMM.summary.proportions.csv", sep = "."))
write.csv(x = s_anns, file = paste(name, "sample.classifications.GMM.csv", sep = "."))


# make density plots of normalized hto counts with sample classifiers
par(mfrow = c(3,2))
for (j in 1:nrow(norm_counts)) {
  plot(density(norm_counts[j,]),
       main = paste(rownames(norm_counts)[j], "HTO CLR Normalized Read Counts", sep = " "))
  # add the mean value for each component
  abline(v = params[[j]]$mean[1], lty = 2, col = "red")
  text(x = params[[j]]$mean[1]+0.5, y = 0.5, labels = paste("Mean: ", round(params[[j]]$mean[1], digits = 2), sep = ""))
  abline(v = params[[j]]$mean[2], lty = 2, col = "red")
  text(x = params[[j]]$mean[2]+0.5, y = 0.75, labels = paste("Mean: ", round(params[[j]]$mean[2], digits = 2), sep = ""))
  
}

