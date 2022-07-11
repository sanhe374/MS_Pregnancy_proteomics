# Correlations of P4 and proteins levels-----------------------------------


# Author: Sandra Hellberg
# email: sandra.hellberg@liu.se

# Date: 220325

# 1. Importation and formatting the input data
# 2. Correlation using cor
# 3. Visualize of significant proteins per time point

# Load packages-----------------------------------------------------------------

pack_R <- c("stats", "dplyr", "pheatmap", "stringr", "corrplot", "ggplot2", "ppcor", "hrbrthemes",
            "readxl", "gridExtra", "ggpubr")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(23)


# 1. Importation and formatting the input data----------------------------------
# meta data
md <- read_excel("Md_Olink.xlsx") # read in meta data file
md <- data.frame(md) # convert to a data frame

# Extract data frame with sample ID and time point
md_timepoint <- md[,c("Sample_ID","TimePoint")]
rownames(md_timepoint) <- md_timepoint$Sample_ID
md_timepoint <- md_timepoint[,-1,drop=FALSE]
md_timepoint <- md_timepoint[order(rownames(md_timepoint)),,drop=FALSE]


# Extract data frame with sample ID and status
md_status <- md[,c("Sample_ID","Individual")]
rownames(md_status) <- md_status$Sample_ID
md_status <- md_status[,-1,drop=FALSE]
md_status <- md_status[order(rownames(md_status)),,drop=FALSE]

# proteins
proteins <- read_excel("proteins_Olink.xlsx") # read in protein data
proteins <- data.frame(proteins)
rownames(proteins) <- proteins$mediator_name
proteins <- proteins[,-1]
proteins <- t(proteins)
proteins <- proteins[order(rownames(proteins)),,drop=FALSE]
proteins <- data.frame(proteins)
proteins_numeric <- sapply(proteins, as.numeric) # convert to numeric
rownames(proteins_numeric) <- rownames(proteins)
proteins <- data.frame(proteins_numeric)
rm(proteins_numeric)

# Progesterone
P4 <- proteins["P4"]
all(rownames(P4) == rownames(proteins))

# Remove P4 and E2 from the protein data
col_names_df_to_remove<-c("P4","E2")
proteins <- proteins[,!(colnames(proteins) %in% col_names_df_to_remove)]


# 2. Correlation using cor------------------------------------------------------

corr_P4 <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(P4$P4, proteins[,i])
  corr_P4[1,i] <- corr[[4]]
  corr_P4[2,i] <- corr[[3]]
}

rownames(corr_P4) <- c("p", "p_val")
colnames(corr_P4) <- colnames(proteins)

corr_P4 <- t(corr_P4)
corr_P4 <- data.frame(corr_P4)

corr_sig_P4 <- corr_P4[abs(corr_P4$r) > 0.5,]

# 3. Visualize of significant proteins ---------------------------
# Visualization of proteins that were differentially expressed and had a 
# significant correlation with P4

proteins_for_graph <- data.frame(proteins=rownames(corr_sig_P4))
p <- list()

for (i in 1:nrow(proteins_for_graph)) {
  prot <- data.frame(protein=proteins[,proteins_for_graph[i,], drop=FALSE])
  all(rownames(P4)==rownames(prot))
  P4_prot <- data.frame(P4 = P4$P4,protein = prot[,1])
  rownames(P4_prot) <- rownames(P4)
  all(rownames(P4_prot)==rownames(md_timepoint))
  P4_prot$TimePoint <- md_timepoint$TimePoint
  
  p[[i]] <- ggplot(P4_prot, aes(x=protein, y=P4, color=TimePoint)) + 
    geom_point(size=6, aes(fill=TimePoint), colour="black", pch=21) +
    geom_smooth(method=lm , color="black", se=FALSE) +
    ggtitle(proteins_for_graph[i,]) + 
    theme_classic() +
    theme(legend.position="none") +
    grids(linetype = "dashed")
  
  # ggsave(paste0(proteins_for_graph[i,], ".pdf"))
  
}

grid.arrange(grobs = p, layout_matrix = rbind(c(1:4), c(5:8), c(9:12),
                                              c(13:16)))

# save plots

# End of script-----------------------------------------------------------------



