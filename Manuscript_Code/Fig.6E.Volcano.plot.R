rm(list = ls())
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(Vennerable)
# load in marker genes
lps_markers <- read.csv("../../Monocyte.doublets.removed.LPS.vs.Unactivated.marker.genes.csv",
                        header = T, check.names = F, stringsAsFactors = F, row.names = 1)

# volcano plots of comparisons (ggrepel to label)
lps_markers$log_pvalue <- -log10(lps_markers$p_val_adj)
id_inf <- which(is.infinite(lps_markers$log_pvalue))
lps_markers$log_pvalue[id_inf] <- 310
# which are significant
id_sig <- which(lps_markers$p_val_adj < 0.05 & abs(lps_markers$avg_logFC) > 0.25)
lps_markers$Significant <- "No"
lps_markers$Significant[id_sig] <- "Yes"

# add labels (top 10 up and down)
lps_markers$Label <- ""
u1 <- strsplit(x = rownames(lps_markers), split = "-")
u2 <- unlist(lapply(u1, `[[`, 1))
id_cite <- which(u2 %in% "CITE")
c1 <- unlist(lapply(u1[id_cite], `[[`, 2))
u2[id_cite] <- paste(u2[id_cite], c1, sep = "_")
u2[1:10]
# match up genes in df
lps_markers$Symbol <- u2
# order matrix
lps_markers_ord <- lps_markers[order(lps_markers$avg_logFC, decreasing = T), ]
lps_markers_ord$Label[1:10] <- lps_markers_ord$Symbol[1:10]
# order again
lps_markers_ord_2 <- lps_markers_ord[order(lps_markers_ord$avg_logFC, decreasing = F), ]
lps_markers_ord_2$Label[1:10] <- lps_markers_ord_2$Symbol[1:10]

p<-ggplot(lps_markers_ord_2, aes(x = avg_logFC, y = log_pvalue, color = Significant, label = Label)) + geom_point(size = 3, alpha =0.5) +
  scale_x_continuous(breaks = seq(from = -1, to = 3, by = 0.5), limits = c(-1, 3)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 310, by = 50), limits = c(0, 310)) +
  scale_color_manual(values = c("grey", "red")) + 
  geom_label_repel(min.segment.length = 0) + 
  labs(x = "Average Log2FC (Stimulated vs. Unactivated)", y = "-log10 FDR adjusted P-value",
       title = "LPS Activated vs. Unactivated Monocytes")
plot(p)

p<-ggplot(lps_markers_ord_2, aes(x = avg_logFC, y = log_pvalue, color = Significant, label = Label)) + geom_point(size = 3, alpha =0.5) +
  scale_x_continuous(breaks = seq(from = -1, to = 3, by = 0.5), limits = c(-1, 3)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 310, by = 50), limits = c(0, 310)) +
  scale_color_manual(values = c("grey", "red")) + 
  labs(x = "Average Log2FC (Stimulated vs. Unactivated)", y = "-log10 FDR adjusted P-value",
       title = "LPS Activated vs. Unactivated Monocytes")
plot(p)

