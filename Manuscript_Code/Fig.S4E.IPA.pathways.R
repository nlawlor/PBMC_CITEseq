library(pheatmap)
library(ggplot2)
library(RColorBrewer)

# read in ipa disease enrichment results (values are already -log10 p-value)
ipa_dis <- read.delim("../../Monocyte.state.genes.comparison.diseases.IPA.txt",
                      header = T, check.names = F, stringsAsFactors = F, skip = 2,
                      row.names = 1)
# remove entries with zero
ipa_dis_sel <- ipa_dis[rowMeans(ipa_dis) > -log10(0.05) ,]
dim(ipa_dis_sel)
colnames(ipa_dis_sel) <- c("Baseline", "LPS2", "LPS3")

# scaling done in pheatmap
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

dis_scale <- ipa_dis_sel
lps2 <- dis_scale[order(dis_scale$LPS2, decreasing = T), ]
lps3 <- dis_scale[order(dis_scale$LPS3, decreasing = T), ]
basel <- dis_scale[order(dis_scale$Baseline, decreasing = T), ]

# make a heatmap of the results?
groups_use <- Reduce(union, list(rownames(lps2)[1:15], rownames(lps3)[1:15],
                                 rownames(basel)[1:15]))
length(groups_use)

mat_res <- ipa_dis_sel[groups_use, ]
# heatmap with significance labels
test_labels <- mat_res
test_labels[test_labels < -log10(0.05)] <- ""
test_labels[test_labels >= -log10(0.05)] <- "*"

breaksList = seq(-15, 15, by = 1)
mat_res$Baseline <- -mat_res$Baseline
ipa_dis_sel$Baseline <- -as.numeric(ipa_dis_sel$Baseline)

# pheatmap of results; Fig. S3E
pheatmap(mat_res, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(length(breaksList)),
         scale = "none", main = "IPA Diseases", display_numbers = test_labels,
         fontsize_number = 12, breaks = breaksList)


