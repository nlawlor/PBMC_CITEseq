library(RColorBrewer)
library(ggplot2)
library(plyr)
library(scales)

# read in donor information
donors <- read.csv("../../CZI.Donor.metadata.csv", header = T, check.names = F, stringsAsFactors = F, row.names = 1)

# flow cytometry data
flow_data <- read.csv("../../Donor.cell.percentages.csv", header = T, check.names = F, stringsAsFactors = F)

# order by females and males
Gender <- rep("Male", nrow(flow_data))
id_females <- which(flow_data$Donor %in% c("11", "2", "6", "7", "10"))
Gender[id_females] <- "Female"

# order plots by females and then males
flow_data$Gender <- Gender

flow_data <- within(flow_data,
                    CellType <- factor(CellType,
                                       levels=c("Bcell", "Monocytes", "Tcell" ,"NK")))


# read in single cell metadata 
dem_meta_single <- readRDS("../../CZI_Fourth_Run_doublets_removed_no_IgM_IgG_metadata.Rds")
dim(dem_meta_single)
sp1 <- strsplit(x = dem_meta_single$BEST.GUESS, split = ",")
dem_meta_single$SNG.POSTERIOR <- unlist(lapply(sp1, `[[`, 1))

# make df of cell numbers per donor
cell_types <- list(c("B"), c("CD14_Mono"), c("CD4T_Naive", "CD4T_Mem", "CD8T_Naive", "CD8T_Mem"), c("NK"))
cell_t_nam <- c("Bcell", "Monocytes", "Tcell", "NK")
dons <- unique(dem_meta_single$SNG.POSTERIOR)
cell_num <- NULL
cell_name <- NULL
don_name <- NULL
pers <- NULL

for (d in 1:length(cell_types)) {
  for (q in 1:length(dons)) {
    df_sel <- dem_meta_single[dem_meta_single$Treatment == "Control" & dem_meta_single$SNG.POSTERIOR == dons[q] & dem_meta_single$Annotation %in% unlist(cell_types[d]), ]
    df_don <- dem_meta_single[dem_meta_single$Treatment == "Control" & dem_meta_single$SNG.POSTERIOR == dons[q] & dem_meta_single$Annotation %in% unlist(cell_types), ]
    cell_num <- c(cell_num, nrow(df_sel))
    cell_name <- c(cell_name, cell_t_nam[d])
    don_name <- c(don_name, dons[q])
    pers <- c(pers, round(100 * (nrow(df_sel)/nrow(df_don))))
  }
}
cell_df <- data.frame(Cell_Number = cell_num, Cell_Type = cell_name, Donor = don_name, Percent = pers)

cell_data <- ddply(cell_df, .(Cell_Number),
                   transform, pos = cumsum(Cell_Number) - (0.5 * Cell_Number))

# reorder levels
cell_data <- within(cell_data,
                    Cell_Type <- factor(Cell_Type,
                                    levels=c("Bcell", "Monocytes", "Tcell", "NK")))

# add donor number to data table
Donor_Number <- as.character(cell_data$Donor)
# match up donor name and number
d1 <- strsplit(x = rownames(donors), split = "10141-")
don_num <- unlist(lapply(d1, `[[`, 2))
donors$Donor_Number <- as.character(as.numeric(don_num))
for (j in 1:nrow(donors)) {
  idx <- which(Donor_Number == donors$Demux_Name[j])
  Donor_Number[idx] <- donors$Donor_Number[j]
}
cell_data$Donor_Number <- Donor_Number

# combine data
colnames(cell_data)[2] <- "CellType"
colnames(cell_data)[4] <- "Control"
colnames(flow_data)[1] <- "Donor_Number"
cell_data$Tech <- "RNA"
flow_data$Tech <- "Flow"
comb_df <- rbind(flow_data[, c(1, 3,4, 6)], cell_data[, c(6, 4, 2, 7)])
comb_df[1:4,]


# make a single X/Y plot of single cell vs. Flow proportions (4 different lines for each cell type)
comb_df_sing <- comb_df[comb_df$Tech == "RNA",]
comb_df_flow <- comb_df[comb_df$Tech == "Flow",]
# need to order matrices by same donor and cell type
id_ord <- NULL
for (o in 1:nrow(comb_df_flow)) {
  idx <- which(comb_df_sing$Donor_Number == comb_df_flow$Donor_Number[o] & comb_df_sing$CellType == comb_df_flow$CellType[o])
  id_ord <- c(id_ord, idx)
}
df_sing_ord <- comb_df_sing[id_ord,]
df_sing_ord[1:4,]
comb_df_flow[1:4,]
comb_df_flow$Single_Cell_Prop <- df_sing_ord$Control
# add gender
Gender <- rep("Male", nrow(comb_df_flow))
id_females <- which(comb_df_flow$Donor_Number %in% c("11", "2", "6", "7", "10"))
Gender[id_females] <- "Female"
comb_df_flow$Sex <- Gender

# plot of cell proportion correlations
p1 <- ggplot(comb_df_flow, aes(y=Control, x=Single_Cell_Prop, group = CellType, color = CellType)) + 
  geom_point(aes(color=CellType, pch = Sex, size = 1.5, alpha = 0.5)) + 
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  labs(y = "Flow Cytometry Cell %", x = "scRNA-seq Cell %") +
  scale_color_manual(values = c("#fdbf6f", "red2", "#2171b5", "#c2a5cf"))
plot(p1)

# pearson R
cor(comb_df_flow$Single_Cell_Prop, comb_df_flow$Control)

