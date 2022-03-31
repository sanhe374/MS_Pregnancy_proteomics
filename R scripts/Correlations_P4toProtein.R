# Correlations of P4 and proteins levels-----------------------------------


# Author: Sandra Hellberg
# email: sandra.hellberg@liu.se

# Date: 220325

# 1. Importation and formatting the input data
# 2. Correlation using cor
# 3. Visualize of significant proteins per time point

# Load packages-----------------------------------------------------------------

pack_R <- c("stats", "dplyr", "pheatmap", "stringr", "corrplot", "ggplot2", "ppcor", "hrbrthemes",
            "readxl")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(23)


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

# Progesterone
P4 <- read.delim("P4_levels.txt", sep = " ")
all(rownames(P4) == rownames(proteins))


# 2. Correlation using cor------------------------------------------------------

# 1st trimester
idx_1st <- grep("_1", rownames(proteins))
proteins_1st <- proteins[idx_1st,]
P4_1st <- P4[idx_1st,,drop=FALSE]
all(rownames(proteins_1st)==rownames(P4_1st))

corr_P4_1st <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(P4_1st$P4, proteins_1st[,i])
  corr_P4_1st[1,i] <- corr[[4]]
  corr_P4_1st[2,i] <- corr[[3]]
}

rownames(corr_P4_1st) <- c("r", "p_val")
colnames(corr_P4_1st) <- colnames(proteins)

corr_P4_1st <- t(corr_P4_1st)
corr_P4_1st <- data.frame(corr_P4_1st)
corr_P4_1st$padjust <- p.adjust(corr_P4_1st$p_val)


corr_P4_1st_sig <- corr_P4_1st[corr_P4_1st$p_val < 0.05,]


# 2nd trimester
idx_2nd <- grep("_2", rownames(proteins))
proteins_2nd <- proteins[idx_2nd,]
P4_2nd <- P4[idx_2nd,,drop=FALSE]
all(rownames(proteins_2nd)==rownames(P4_2nd))

corr_P4_2nd <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(P4_2nd$P4, proteins_2nd[,i])
  corr_P4_2nd[1,i] <- corr[[4]]
  corr_P4_2nd[2,i] <- corr[[3]]
}

rownames(corr_P4_2nd) <- c("r", "p_val")
colnames(corr_P4_2nd) <- colnames(proteins)

corr_P4_2nd <- t(corr_P4_2nd)
corr_P4_2nd <- data.frame(corr_P4_2nd)

corr_P4_2nd_sig <- corr_P4_2nd[corr_P4_2nd$p_val < 0.05,]

# 3rd trimester
idx_3rd <- grep("_3", rownames(proteins))
proteins_3rd <- proteins[idx_3rd,]
P4_3rd <- P4[idx_3rd,,drop=FALSE]
all(rownames(proteins_3rd)==rownames(P4_3rd))

corr_P4_3rd <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(P4_3rd$P4, proteins_3rd[,i])
  corr_P4_3rd[1,i] <- corr[[4]]
  corr_P4_3rd[2,i] <- corr[[3]]
}

rownames(corr_P4_3rd) <- c("r", "p_val")
colnames(corr_P4_3rd) <- colnames(proteins)

corr_P4_3rd <- t(corr_P4_3rd)
corr_P4_3rd <- data.frame(corr_P4_3rd)

corr_P4_3rd_sig <- corr_P4_3rd[corr_P4_3rd$p_val < 0.05,]

# PP 
idx_PP <- grep("_4", rownames(proteins))
proteins_PP <- proteins[idx_PP,]
P4_PP <- P4[idx_PP,,drop=FALSE]
all(rownames(proteins_PP)==rownames(P4_PP))

corr_P4_PP <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(P4_PP$P4, proteins_PP[,i])
  corr_P4_PP[1,i] <- corr[[4]]
  corr_P4_PP[2,i] <- corr[[3]]
}

rownames(corr_P4_PP) <- c("r", "p_val")
colnames(corr_P4_PP) <- colnames(proteins)

corr_P4_PP <- t(corr_P4_PP)
corr_P4_PP <- data.frame(corr_P4_PP)

corr_P4_PP_sig <- corr_P4_PP[corr_P4_PP$p_val < 0.05,]


# 3. Visualize of significant proteins per time point---------------------------
# Visualization of proteins that were differentially expressed and had a 
# significant correlation with P4

# LIFR 3rd trimester P4
proteins_LIFR <- data.frame(LIF_R=proteins[,"LIF_R", drop=FALSE])
timepoint_3rd <- grep("_3", rownames(proteins))
LIFR_3rd <- proteins_LIFR[timepoint_3rd,,drop=FALSE]
P4_3rd <- P4[timepoint_3rd,,drop=FALSE]
md_3rd <- md_status[timepoint_3rd,,drop=FALSE]

all(rownames(P4_3rd)==rownames(LIFR_3rd))

P4_LIFR <- data.frame(P4 = P4_3rd$P4,LIFR = LIFR_3rd$LIF_R)

ggplot(P4_LIFR, aes(x=LIFR, y=P4, color=md_3rd$Individual)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)+
  theme_ipsum()


# CDCP1 3rd trimester P4
proteins_CDCP1 <- data.frame(CDCP1=proteins[,"CDCP1", drop=FALSE])
CDCP1_3rd <- proteins_CDCP1[timepoint_3rd,,drop=FALSE]

all(rownames(P4_3rd)==rownames(CDCP1_3rd))

P4_CDCP1 <- data.frame(P4 = P4_3rd$P4,CDCP1 = CDCP1_3rd$CDCP1)

ggplot(P4_CDCP1, aes(x=CDCP1, y=P4, color=md_3rd$Individual)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)+
  theme_ipsum()


# HGF 3rd trimester P4
proteins_HGF <- data.frame(HGF=proteins[,"HGF", drop=FALSE])
HGF_3rd <- proteins_HGF[timepoint_3rd,,drop=FALSE]

all(rownames(P4_3rd)==rownames(HGF_3rd))

P4_HGF <- data.frame(P4 = P4_3rd$P4,HGF = HGF_3rd$HGF)

ggplot(P4_HGF, aes(x=HGF, y=P4, color=md_3rd$Individual)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)+
  theme_ipsum()


# IL18R1 3rd trimester P4
proteins_IL18R1 <- data.frame(IL18R1=proteins[,"IL_18R1", drop=FALSE])
IL18R1_3rd <- proteins_IL18R1[timepoint_3rd,,drop=FALSE]

all(rownames(P4_3rd)==rownames(IL18R1_3rd))

P4_IL18R1 <- data.frame(P4 = P4_3rd$P4,IL18R1= IL18R1_3rd $IL_18R1)

ggplot(P4_IL18R1, aes(x=IL18R1, y=P4, color=md_3rd$Individual)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)+
  theme_ipsum()


# End of script-----------------------------------------------------------------



