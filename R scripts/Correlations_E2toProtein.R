# Correlations of E2 and proteins levels-----------------------------------


# Author: Sandra Hellberg
# email: sandra.hellberg@liu.se

# Date: 220325

# 1. Importation and formatting the input data
# 2. Formatting and subsetting input data
# 3. Visualize of significant proteins per time point

# Load packages-----------------------------------------------------------------

pack_R <- c("stats", "dplyr", "pheatmap", "stringr", "corrplot", "ggplot2", "ppcor", "hrbrthemes",
            "readxl")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(150)


# 1. Importation and formatting the input data----------------------------------
# meta data
md <- read_excel("Md_Olink.xlsx") # read in meta data file
md <- data.frame(md) # convert to a data frame
md$Sample_ID <- gsub("P", "", md$Sample_ID)


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
rownames(proteins) <- gsub("P", "", rownames(proteins)) # remove P before sampleID
proteins <- proteins[order(rownames(proteins)),,drop=FALSE]
proteins <- data.frame(proteins)
proteins_numeric <- sapply(proteins, as.numeric) # convert to numeric
rownames(proteins_numeric) <- rownames(proteins)
proteins <- data.frame(proteins_numeric)
rm(proteins_numeric)

# Estrogen
E2 <- read.delim("E2_levels.txt", sep = " ")
all(rownames(E2) == rownames(proteins))


# 2. Correlation using cor------------------------------------------------------

# 1st trimester
idx_1st <- grep("_1", rownames(proteins))
proteins_1st <- proteins[idx_1st,]
E2_1st <- E2[idx_1st,,drop=FALSE]
all(rownames(proteins_1st)==rownames(E2_1st))

corr_E2_1st <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(E2_1st$E2, proteins_1st[,i])
  corr_E2_1st[1,i] <- corr[[4]]
  corr_E2_1st[2,i] <- corr[[3]]
}

rownames(corr_E2_1st) <- c("r", "p_val")
colnames(corr_E2_1st) <- colnames(proteins)

corr_E2_1st <- t(corr_E2_1st)
corr_E2_1st <- data.frame(corr_E2_1st)
corr_E2_1st$padjust <- p.adjust(corr_E2_1st$p_val)


corr_E2_1st_sig <- corr_E2_1st[corr_E2_1st$p_val < 0.05,]


# 2nd trimester
idx_2nd <- grep("_2", rownames(proteins))
proteins_2nd <- proteins[idx_2nd,]
E2_2nd <- E2[idx_2nd,,drop=FALSE]
all(rownames(proteins_2nd)==rownames(E2_2nd))

corr_E2_2nd <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(E2_2nd$E2, proteins_2nd[,i])
  corr_E2_2nd[1,i] <- corr[[4]]
  corr_E2_2nd[2,i] <- corr[[3]]
}

rownames(corr_E2_2nd) <- c("r", "p_val")
colnames(corr_E2_2nd) <- colnames(proteins)

corr_E2_2nd <- t(corr_E2_2nd)
corr_E2_2nd <- data.frame(corr_E2_2nd)
corr_E2_2nd$padjust <- p.adjust(corr_E2_2nd$p_val)


corr_E2_2nd_sig <- corr_E2_2nd[corr_E2_2nd$p_val < 0.05,]

# 3rd trimester
idx_3rd <- grep("_3", rownames(proteins))
proteins_3rd <- proteins[idx_3rd,]
E2_3rd <- E2[idx_3rd,,drop=FALSE]
all(rownames(proteins_3rd)==rownames(E2_3rd))

corr_E2_3rd <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(E2_3rd$E2, proteins_3rd[,i])
  corr_E2_3rd[1,i] <- corr[[4]]
  corr_E2_3rd[2,i] <- corr[[3]]
}

rownames(corr_E2_3rd) <- c("r", "p_val")
colnames(corr_E2_3rd) <- colnames(proteins)

corr_E2_3rd <- t(corr_E2_3rd)
corr_E2_3rd <- data.frame(corr_E2_3rd)
corr_E2_3rd$padjust <- p.adjust(corr_E2_3rd$p_val)



corr_E2_3rd_sig <- corr_E2_3rd[corr_E2_3rd$p_val < 0.05,]

# PP 
idx_PP <- grep("_4", rownames(proteins))
proteins_PP <- proteins[idx_PP,]
E2_PP <- E2[idx_PP,,drop=FALSE]
all(rownames(proteins_PP)==rownames(E2_PP))

corr_E2_PP <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(E2_PP$E2, proteins_PP[,i])
  corr_E2_PP[1,i] <- corr[[4]]
  corr_E2_PP[2,i] <- corr[[3]]
}

rownames(corr_E2_PP) <- c("r", "p_val")
colnames(corr_E2_PP) <- colnames(proteins)

corr_E2_PP <- t(corr_E2_PP)
corr_E2_PP <- data.frame(corr_E2_PP)
corr_E2_PP$padjust <- p.adjust(corr_E2_PP$p_val)

corr_E2_PP_sig <- corr_E2_PP[corr_E2_PP$p_val < 0.05,]


# 3. Visualize of significant proteins per time point---------------------------
# Visualization of proteins that were differentially expressed and had a 
# significant correlation with E2

# CSF1 3rd trimester E2
proteins_CSF1 <- data.frame(CSF1=proteins[,"CSF_1", drop=FALSE])
timepoint_3rd <- grep("_3", rownames(proteins))
CSF1_3rd <- proteins_CSF1[timepoint_3rd,,drop=FALSE]
E2_3rd <- E2[timepoint_3rd,,drop=FALSE]
md_3rd <- md_status[timepoint_3rd,,drop=FALSE]


all(rownames(E2_3rd)==rownames(CSF1_3rd))

E2_CSF1 <- data.frame(E2 = E2_3rd$E2,CSF1 = CSF1_3rd$CSF_1)

ggplot(E2_CSF1, aes(x=CSF1, y=E2,color=md_3rd$Individual)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()

# End of script-----------------------------------------------------------------