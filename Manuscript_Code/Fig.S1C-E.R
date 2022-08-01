# Compare demuxlet doublets and HTO doublets
rm(list = ls())
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
# demuxlet results
best_all <- readRDS("../../CZI.Fourth.Run.all.barcodes.Popscle.Demuxlet.BEST.Rds")
dim(best_all)
# hto data
gmm <- read.csv("../../CZI.PMBC.HTO.pooled.sample.classifications.GMM.csv",
                header = T, check.names = F, stringsAsFactors = F, row.names = 2)
dim(gmm)
rownames(gmm) <- gsub(x = rownames(gmm), pattern = "HTO", replacement = "Sample")
length(which(rownames(gmm) %in% rownames(best_all)))

# remove anti-igM/IgG stim cells
id_igm <- which(grepl(x = gmm$Sample, pattern = "Bcell"))
length(id_igm)
gmm <- gmm[-id_igm, ]
dim(gmm)

int_names <- intersect(rownames(best_all), rownames(gmm))
dem_sel <- best_all[int_names, ]
gmm_sel <- gmm[int_names, ]
table(rownames(dem_sel) == rownames(gmm_sel))

# add metadata to dem columns
dem_ext <- cbind(dem_sel, gmm_sel)

tab_hto <- dem_ext$Sample
id_quad <- which(grepl(x = tab_hto, glob2rx("*\\+*\\+*\\+*")))
id_tr1 <- which(grepl(x = tab_hto, glob2rx("*\\+*\\+*")))
id_trip <- setdiff(id_tr1, id_quad)
id_d1 <- which(grepl(x = tab_hto, glob2rx("*\\+*")))
id_doub <- setdiff(id_d1, c(id_quad, id_trip))
table(tab_hto[id_doub])

# vector of hto type
HTO_Type <- rep("Singlet", length(tab_hto))
HTO_Type[id_doub] <- "Doublet"
HTO_Type[id_trip] <- "Triplet"
HTO_Type[id_quad] <- "Quadruplet"
table(HTO_Type)
id_no <- which(tab_hto == "None")
HTO_Type[id_no] <- "None"
table(HTO_Type)
dem_ext$HTO_Droplet_Type <- HTO_Type

hto_table <- table(dem_ext$Sample)
hto_df <- data.frame(HTO_Type = names(hto_table), Number = as.numeric(hto_table))
Droplet_Type <- rep("Singlet", nrow(hto_df))
id_q <- which(grepl(x = hto_df$HTO_Type, glob2rx("*\\+*\\+*\\+*")))
id_t <- which(grepl(x = hto_df$HTO_Type, glob2rx("*\\+*\\+*")))
id_tr <- setdiff(id_t, id_q)
id_d1 <- which(grepl(x = hto_df$HTO_Type, glob2rx("*\\+*")))
id_db <- setdiff(id_d1, c(id_q, id_t))
Droplet_Type[id_q] <- "Quadruplet"
Droplet_Type[id_tr] <- "Triplet"
Droplet_Type[id_db] <- "Doublet"
id_none <- which(hto_df$HTO_Type == "None")
Droplet_Type[id_none] <- "Empty"
hto_df$Droplet_Type <- Droplet_Type


# plot cell numbers of singlets, multiplets, and empty (Fig. S1C)
hto_ident <- dem_ext$Sample_Short
id_sing <- which(hto_ident %in% c("Control", "Mono", "Tcell"))
hto_ident[id_sing] <- "Singlet"
hto_ident <- gsub(hto_ident, pattern = "None", replacement = "Empty")
hto_ident_df <- as.data.frame(table(hto_ident))
p1 <- ggplot(hto_ident_df, aes(x=reorder(hto_ident, -Freq), y=Freq)) + 
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(p1)

# make a plot of cell identities (by demuxlet) after removing doublets from cell hashing no igM igG
dem_no_hto <- dem_ext[dem_ext$HTO_Droplet_Type == "Singlet", ]
dim(dem_no_hto)

dem_no_hto_tab <- as.data.frame(table(dem_no_hto$DROPLET.TYPE))
# Fig S1D
p1 <- ggplot(dem_no_hto_tab, aes(x=reorder(Var1, -Freq), y=Freq)) + 
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Demuxlet Identity", y = "Number of Cells")
plot(p1)


# vector of HTO identity
HTO_ident <- dem_ext$HTO_Droplet_Type
for (i in 1:length(HTO_ident)) {
  # if not singlet, skip
  if (HTO_ident[i] == "Singlet") {
    HTO_ident[i] <- as.character(dem_ext$Sample_Short[i])
  } else {}
}
dem_ext$HTO_ident <- HTO_ident

dem_hto <- as.data.frame(table(dem_ext$HTO_ident, dem_ext$DROPLET.TYPE))
colnames(dem_hto) <- c("HTO_Ident", "Demuxlet_Ident", "Number")
sum(dem_hto$Number)
dim(dem_ext)

# make a plot by donor and HTO identity (6 bars by HTO identity (3 treatments, empty, doublet, triplet), split by demuxlet)
dem_hto <- within(dem_hto,
                  HTO_Ident <- factor(HTO_Ident,
                                  levels=c("Control", "Mono", "Tcell", "None", "Doublet", "Triplet")))

# Fig S1E
p1 <- ggplot(dem_hto, aes(fill=Demuxlet_Ident, x=HTO_Ident, y = Number)) + 
  geom_bar(stat="identity", position="fill") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = "Cell Proportions") + scale_fill_manual(values = c("#a2a637", "#f3766e", "#38ace2"))
plot(p1)

