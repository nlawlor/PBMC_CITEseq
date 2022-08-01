# Calculate cell type proportions for each
rm(list = ls())
library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse)

# load in CD8T subpop info
cd8 <- readRDS("../../CD8T.doublets.removed.no.IgM.seurat.Rds")
cd8_meta <- cd8@meta.data
cd8_meta <- cd8_meta[cd8_meta$Treatment == "Baseline",]
dim(cd8_meta)
tab_don <- as.data.frame(table(cd8_meta$New_Annotations, cd8_meta$BEST.GUESS))
tab_don[1:4,]

treat_cols <- c("#F3766E", "#7CAF42", "#2278BE", "#A781BA")
gp1 <- ggplot(tab_don, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") + labs(x = "", y = "# of cells per Donor") +
  scale_fill_manual(values = treat_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot(gp1)

gp1 <- ggplot(tab_don, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", position = "fill") + labs(x = "", y = "% of cells per Donor") +
  scale_fill_manual(values = treat_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot(gp1)


# load in CD8T subpop info
cd4 <- readRDS("../../CD4T.doublets.removed.no.IgM.seurat.Rds")
cd4_meta <- cd4@meta.data
cd4_meta <- cd4_meta[cd4_meta$Treatment == "Baseline",]
cd4_meta <- cd4_meta[cd4_meta$New_Annotations != "CD4T_Naive_Actv", ]
dim(cd4_meta)
tab_don <- as.data.frame(table(cd4_meta$New_Annotations, cd4_meta$BEST.GUESS))
tab_don[1:4,]

treat_cols <- c("#F3766E", "#BFBF41", "#7CAF42", "#38ABE2", "#A781BA")
gp1 <- ggplot(tab_don, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") + labs(x = "", y = "# of cells per Donor") +
  scale_fill_manual(values = treat_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot(gp1)

gp1 <- ggplot(tab_don, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", position = "fill") + labs(x = "", y = "% of cells per Donor") +
  scale_fill_manual(values = treat_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot(gp1)

