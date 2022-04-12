# Correlations of E2 and proteins levels-----------------------------------


# Author: Sandra Hellberg
# email: sandra.hellberg@liu.se

# Date: 220325

# 1. Importation and formatting the input data
# 2. Formatting and subsetting input data
# 3. Visualize of significant proteins per time point

# Load packages-----------------------------------------------------------------

pack_R <- c("stats", "dplyr", "pheatmap", "stringr", "corrplot", "ggplot2", "ppcor", "hrbrthemes",
            "readxl", "gridExtra", "ggpubr")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(150)


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

# Estrogen
E2 <- proteins["E2"]
all(rownames(E2) == rownames(proteins))

# Remove P4 and E2 from the protein data
col_names_df_to_remove<-c("P4","E2")
proteins <- proteins[,!(colnames(proteins) %in% col_names_df_to_remove)]

# 2. Correlation using cor------------------------------------------------------

corr_E2 <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(E2$E2, proteins[,i])
  corr_E2[1,i] <- corr[[4]]
  corr_E2[2,i] <- corr[[3]]
}

rownames(corr_E2) <- c("r", "p_val")
colnames(corr_E2) <- colnames(proteins)

corr_E2 <- t(corr_E2)
corr_E2 <- data.frame(corr_E2)

corr_sig_E2 <- corr_E2[abs(corr_E2$r) > 0.5,]

# 3. Visualize of significant proteins per time point---------------------------
# Visualization of proteins that were differentially expressed and had a 
# significant correlation with E2

proteins_for_graph <- data.frame(proteins=rownames(corr_sig_E2))
p <- list()

for (i in 1:nrow(proteins_for_graph)) {
  prot <- data.frame(protein=proteins[,proteins_for_graph[i,], drop=FALSE])
  all(rownames(E2)==rownames(prot))
  E2_prot <- data.frame(E2 = E2$E2,protein = prot[,1])
  rownames(E2_prot) <- rownames(E2)
  all(rownames(E2_prot)==rownames(md_timepoint))
  E2_prot$TimePoint <- md_timepoint$TimePoint
  
  p[[i]] <- ggplot(E2_prot, aes(x=protein, y=E2, color=TimePoint)) + 
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

# save p

# End of script-----------------------------------------------------------------