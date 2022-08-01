# Calculate cell type proportions for each
rm(list = ls())
library(ggplot2)
library(ggrepel)
# read in metadata file
meta_comb <- readRDS("../../CZI_Fourth_Run_doublets_removed_no_IgM_IgG_metadata.Rds")

# colors
cols <- c("#fdbf6f", "#ff7f00", "#f7fcb9", "#addd8e", "#41ab5d", "#006837",
          "lightskyblue1", "lightskyblue3", "#2171b5", "#08306b", 
          "rosybrown1", "red2", "red4", "#c2a5cf", "#762a83")
celltype <- c("B", "Actv_B", "CD4T_Naive", "Actv_CD4T_Naive", "CD4T_Mem", "Actv_CD4T_Mem",
              "CD8T_Naive", "Actv_CD8T_Naive", "CD8T_Mem", "Actv_CD8T_Mem",
              "CD16_Mono", "CD14_Mono", "Actv_CD14_Mono", "NK", "Actv_NK")

# cell proportions per condition
meta_sel <- meta_comb[meta_comb$Treatment == "Control", ]
df_treat_cell <- as.data.frame(table(meta_sel$Annotation, meta_sel$BEST.GUESS))
df_treat_cell <- within(df_treat_cell,
                        Var1 <- factor(Var1,
                                      levels=c("B", "CD4T_Naive","CD4T_Mem",
                                               "CD8T_Naive", "CD8T_Mem",
                                               "CD14_Mono", "NK")))
id_c <- which(celltype %in% df_treat_cell$Var1)

# cell numbers at baseline
gp1 <- ggplot(df_treat_cell, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") + labs(x = "", y = "# of cells per Donor", fill = "Cell-type") +
  scale_fill_manual(values = cols[id_c], labels = celltype[id_c]) +
  theme(axis.text.x = element_blank()) + 
  ggtitle("Baseline")
plot(gp1)

# cell percentages at baseline
gp1 <- ggplot(df_treat_cell, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", position = "fill") + labs(x = "", y = "% of cells per Donor", fill = "Cell-type") +
  scale_fill_manual(values = cols[id_c], labels = celltype[id_c]) +
  theme(axis.text.x = element_blank()) +
  ggtitle("Baseline")
plot(gp1)

# cell numbers at anti-CD3/CD28
meta_sel <- meta_comb[meta_comb$Treatment == "Tcell", ]
df_treat_cell <- as.data.frame(table(meta_sel$Annotation, meta_sel$BEST.GUESS))
df_treat_cell <- within(df_treat_cell,
                        Var1 <- factor(Var1,
                                       levels=c("B", "CD4T_Naive","CD4T_Mem",
                                                "CD8T_Naive", "CD8T_Mem",
                                                "CD14_Mono", "NK")))
id_c <- which(celltype %in% df_treat_cell$Var1)

gp1 <- ggplot(df_treat_cell, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") + labs(x = "", y = "# of cells per Donor", fill = "Cell-type") +
  scale_fill_manual(values = cols[id_c], labels = celltype[id_c]) +
  theme(axis.text.x = element_blank()) + 
  ggtitle("anti-CD3/CD28")
plot(gp1)

# cell percentages at baseline
gp1 <- ggplot(df_treat_cell, mapping = aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", position = "fill") + labs(x = "", y = "% of cells per Donor", fill = "Cell-type") +
  scale_fill_manual(values = cols[id_c], labels = celltype[id_c]) +
  theme(axis.text.x = element_blank()) +
  ggtitle("anti-CD3/CD28")
plot(gp1)

