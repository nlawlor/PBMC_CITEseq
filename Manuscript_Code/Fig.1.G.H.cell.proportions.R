# Calculate cell type proportions for each
rm(list = ls())
library(ggplot2)
library(ggrepel)
# read in metadata file
meta_comb <- readRDS("../../CZI_Fourth_Run_doublets_removed_no_IgM_IgG_metadata.Rds")

# colors for treatment conditions
treat_cols <- c("#e41a1c", "#4daf4a", "#984ea3")

# data frane of # of cells by donor and treatment
donor_df <- as.data.frame(table(meta_comb$BEST.GUESS, meta_comb$Treatment))
Updated_Treatment <- donor_df$Var2
Updated_Treatment <- gsub(Updated_Treatment, pattern = "Tcell", replacement = "Anti-CD3/CD28")
Updated_Treatment <- gsub(Updated_Treatment, pattern = "Mono", replacement = "LPS")
Updated_Treatment <- gsub(Updated_Treatment, pattern = "Control", replacement = "Baseline")
donor_df$Treatment <- Updated_Treatment
donor_df[1:4,]

# re-order the levels
donor_df <- within(donor_df,
                   Treatment <- factor(Treatment,
                                      levels=c("Baseline", "LPS", "Anti-CD3/CD28")))
sp1 <- strsplit(x = as.character(donor_df$Var1), split = ",")
donor_df$Donor <- unlist(lapply(sp1, `[[`, 1))

# make plot of cell numbers and proportions
gp1 <- ggplot(donor_df, mapping = aes(x = Donor, y = Freq, fill = Treatment)) + 
  geom_bar(stat = "identity") + labs(x = "", y = "# of cells per Donor") +
  scale_fill_manual(values = treat_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot(gp1)

gp1 <- ggplot(donor_df, mapping = aes(x = Donor, y = Freq, fill = Treatment)) + 
  geom_bar(stat = "identity", position = "fill") + labs(x = "", y = "% of cells per Donor") +
  scale_fill_manual(values = treat_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot(gp1)

