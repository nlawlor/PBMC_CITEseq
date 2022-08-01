# Quantiy proportion of cells that are + for ADT expression
rm(list = ls())
library(mclust)
library(Matrix)
library(pheatmap)
library(Biobase)
library(Seurat)
library(ggplot2)
library(tidyr)
name <- "CZI_Fourth_Run_ADT_doublets_removed"

# read in adt data
adt_counts <- readRDS("../../CZI.PMBC.ADT.pooled.combined.raw.matrix.Rds")
colnames(adt_counts) <- gsub(x = colnames(adt_counts), pattern = "ADT", replacement = "Sample")

# find intersection of all cell barcodes
dem_meta_single <- readRDS("../../CZI_Fourth_Run_doublets_removed_no_IgM_IgG_metadata.Rds")
dim(dem_meta_single)
sp1 <- strsplit(x = dem_meta_single$BEST.GUESS, split = ",")
Donor_Info <- sapply(sp1, `[[`, 1)
dem_meta_single$Donor_Info <- Donor_Info

int_cells <- Reduce(intersect, list(rownames(dem_meta_single), colnames(adt_counts)))
length(int_cells)
adt_sel <- adt_counts[, int_cells]
table(colnames(adt_sel) == rownames(dem_meta_single))

# remove bad struct, total_reads, no_match
ids_adt <- which(rownames(adt_sel) %in% c("bad_struct", "no_match", "total_reads"))
adt_sel <- adt_sel[-ids_adt,]
dim(adt_sel)
adt_sel[1:4,1:4]

# get percentage of B cells that are activated
ser <- CreateSeuratObject(counts = adt_sel)
ser <- NormalizeData(ser, normalization.method = "CLR")
adt_norm <- GetAssayData(ser, slot = "data")

id_cd69 <- which(grepl(x = rownames(adt_norm), pattern = "CD69"))
id_cd25 <- which(grepl(x = rownames(adt_norm), pattern = "CD25"))

# which t cells are activated
id_t <- which(grepl(x = dem_meta_single$Annotation, pattern = "T_") & dem_meta_single$Treatment %in% c("Control", "Tcell"))
cd69 <- Mclust(adt_norm[id_cd69, id_t], G = 2)
cd25 <- Mclust(adt_norm[id_cd25, id_t], G = 2)

# get all treated with cd3 and baseline
meta_cd3_all <- dem_meta_single[dem_meta_single$Treatment %in% c("Tcell"), ]
meta_cont_all <- dem_meta_single[dem_meta_single$Treatment %in% c("Control"), ]
meta_cd3_t <- dem_meta_single[dem_meta_single$Treatment %in% c("Tcell") & grepl(x = dem_meta_single$Annotation, pattern = "T_"), ]
meta_cont_t <- dem_meta_single[dem_meta_single$Treatment %in% c("Control") & grepl(x = dem_meta_single$Annotation, pattern = "T_"), ]
# determine which cells are + for CD69/CD25 in all matrices, then get proportions by donor
length(which(as.numeric(adt_norm[id_cd69, rownames(meta_cd3_all)]) > as.numeric(cd69$parameters$mean[2])))
cd3_all_pos <- which(as.numeric(adt_norm[id_cd69, rownames(meta_cd3_all)]) > as.numeric(cd69$parameters$mean[2]) & as.numeric(adt_norm[id_cd25, rownames(meta_cd3_all)]) > as.numeric(cd25$parameters$mean[2]))
cd3_all_pos_df <- meta_cd3_all[cd3_all_pos, ]
don_cd3_pos <- table(cd3_all_pos_df$Donor_Info)
# table of cells per donor for all cd3
don_cd3_all <- table(meta_cd3_all$Donor_Info)
don_all_cells_cd3 <- round(100* (don_cd3_pos/don_cd3_all), digits = 2)

cont_all_pos <- which(as.numeric(adt_norm[id_cd69, rownames(meta_cont_all)]) > as.numeric(cd69$parameters$mean[2]) & as.numeric(adt_norm[id_cd25, rownames(meta_cont_all)]) > as.numeric(cd25$parameters$mean[2]))
cont_all_pos_df <- meta_cont_all[cont_all_pos, ]
don_cont_pos <- table(cont_all_pos_df$Donor_Info)
# table of cells per donor for all cont
don_cont_all <- table(meta_cont_all$Donor_Info)
don_all_cells_cont <- round(100* (don_cont_pos/don_cont_all), digits = 2)
don_all_cells_cont <- don_cont_all
don_all_cells_cont[1:10] <- 0

# table of cells per donor for b cells at cd3
length(which(as.numeric(adt_norm[id_cd69, rownames(meta_cd3_t)]) > as.numeric(cd69$parameters$mean[2])))
cd3_t_pos <- which(as.numeric(adt_norm[id_cd69, rownames(meta_cd3_t)]) > as.numeric(cd69$parameters$mean[2]) & as.numeric(adt_norm[id_cd25, rownames(meta_cd3_t)]) > as.numeric(cd25$parameters$mean[2]))
cd3_t_pos_df <- meta_cd3_t[cd3_t_pos, ]
don_cd3_t_pos <- table(cd3_t_pos_df$Donor_Info)
# table of b cells per donor for cd3
don_cd3_t <- table(meta_cd3_t$Donor_Info)
don_t_cells_cd3 <- round(100* (don_cd3_t_pos/don_cd3_t), digits = 2)

# table of b cells per donor at control
length(which(as.numeric(adt_norm[id_cd69, rownames(meta_cont_t)]) > as.numeric(cd69$parameters$mean[2])))
cont_t_pos <- which(as.numeric(adt_norm[id_cd69, rownames(meta_cont_t)]) > as.numeric(cd69$parameters$mean[2]) & as.numeric(adt_norm[id_cd25, rownames(meta_cont_t)]) > as.numeric(cd25$parameters$mean[2]))
cont_t_pos_df <- meta_cont_t[cont_t_pos, ]
don_cont_t_pos <- table(cont_t_pos_df$Donor_Info)
# table of b cells per donor for cd3
don_cont_t <- table(meta_cont_t$Donor_Info)
don_t_cells_cont <- round(100* (don_cont_t_pos/don_cont_t), digits = 2)
don_t_cells_cont <- don_cont_t
# values are zero
don_t_cells_cont[1:10] <- 0


boxplot(as.numeric(don_t_cells_cd3), as.numeric(don_t_cells_cont),
        as.numeric(don_all_cells_cd3), as.numeric(don_all_cells_cont), ylim = c(0,100),
        ylab = "% CD69+ in ADT", names = c("T cells (CD3)", "T cells (Control)", "All Cells (CD3)", "All Cells (Control)"),
        main = "")
# make fancy ggplot image
b_cell_df <- data.frame(Percent = c(as.numeric(don_t_cells_cd3), as.numeric(don_t_cells_cont),
                                    as.numeric(don_all_cells_cd3), as.numeric(don_all_cells_cont)),
                        Donor = c(names(don_t_cells_cd3), names(don_t_cells_cont),
                                  names(don_all_cells_cd3), names(don_all_cells_cont)),
                        Group = c(rep("T cells", 20), rep("All cells", 20)),
                        Category = c(rep("T cells (CD3)", 10), rep("T cells (Control)", 10),
                                     rep("All cells (CD3)", 10), rep("All cells (Control)", 10)))

dim(b_cell_df)

b_cell_df <- within(b_cell_df,
                    Category <- factor(Category,
                                    levels=c("T cells (CD3)", "T cells (Control)",
                                             "All cells (CD3)", "All cells (Control)" )))


# make a scatter plot of b cell activation props in flow and ADT
flow_b <- read.csv("../../Tcell.activation.summary.CD69.CD25.csv", header = T, check.names = F, stringsAsFactors = F)
dim(flow_b)
flow_b[1:4,]
flow_b_df <- gather(flow_b, key = "Category", value = "Flow_Percent", `T cells (CD3)`:`All cells (Control)`)
flow_b_df[1:4,]
dim(flow_b_df)
# combine two matrices
both_b <- merge(x = flow_b_df, y = b_cell_df, by.x = c("Donor", "Category"), by.y = c("Donor", "Category"))

p1 <- ggplot(both_b, aes(y=Flow_Percent, x=Percent, group = Category, color = Category)) + 
  geom_point(aes(color=Category, pch = Gender, labels = Donor), size = 4) + 
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(y = "% CD69+/CD25+ in Flow Cytometry", x = "% CD69+/CD25+ in ADT") +
  ggtitle("Flow Cytometry vs. ADT Anti-CD3+CD28 T-cell activation") +
  annotate("text", x=25, y=60, label= paste("R = ", round(cor(both_b$Flow_Percent, both_b$Percent), digits = 2), sep = ""))
plot(p1)

