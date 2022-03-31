# Correlations of P4 and proteins levels-----------------------------------


# Author: Sandra Hellberg
# email: sandra.hellberg@liu.se

# Date: 220325

# 1. Importation and formatting the input data
# 2. Formatting and subsetting input data

# Load packages-----------------------------------------------------------------

pack_R <- c("stats", "dplyr", "pheatmap", "stringr", "corrplot", "ggplot2", "ppcor", "hrbrthemes",
            "readxl")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(23)


# 1. Importation and formatting the input data----------------------------------
md <- read_excel("Md_Olink.xlsx") # read in meta data file
md <- data.frame(md) # convert to a data frame

proteins <- read_excel("proteins_Olink.xlsx") # read in protein data
proteins <- data.frame(proteins)
rownames(proteins) <- proteins$mediator_name # set protein name to row names
proteins <- proteins[,-1]    
proteins_numeric <- sapply(proteins, as.numeric) # convert to numeric
rownames(proteins_numeric) <- rownames(proteins)
proteins <- data.frame(proteins_numeric)

rm(proteins_numeric)


# 2. Formatting and subsetting input data
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

samples_included <- rownames(proteins) # create a vector of all samples included for olink

hormones <- read.delim("Hormone_levels.txt")
colnames(hormones)[1] <- "SampleID"
rownames(hormones) <- hormones$SampleID
hormones <- hormones[,2:3,drop=FALSE]

hormones_olink <- hormones[rownames(hormones) %in% samples_included,]

E2 <- hormones_olink[, 1, drop = F]
E2 <- E2[order(rownames(E2)),,drop=FALSE]

P4 <- hormones_olink[, 2, drop = F]
P4 <- P4[order(rownames(P4)),,drop=FALSE]


all(rownames(P4)==rownames(E2)) # double check that the samples are in the corrected order

# Correlation using all samples-------------------------------------------------

# Correlations with P4
corr_P4 <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  corr <- cor.test(P4$P4, proteins[,i])
  corr_P4[1,i] <- corr[[4]]
  corr_P4[2,i] <- corr[[3]]
}

rownames(corr_P4) <- c("r", "p_val")
colnames(corr_P4) <- colnames(proteins)

corr_P4 <- t(corr_P4)
corr_P4 <- data.frame(corr_P4)

corr_sig_P4 <- corr_P4[abs(corr_P4$r) > 0.5,]

# Correlations with E2

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


# Partial correlation-----------------------------------------------------------

# Progesterone
P4$Timepoint <- md_timepoint$TimePoint
P4$Timepoint[P4$Timepoint=="1st"] <- "1"
P4$Timepoint[P4$Timepoint=="2nd"] <- "2"
P4$Timepoint[P4$Timepoint=="3rd"] <- "3"
P4$Timepoint[P4$Timepoint=="Postpartum"] <- "4"
P4$Timepoint[P4$Timepoint=="Before"] <- "0"
P4$Timepoint[P4$Timepoint=="Non-pregnant"] <- "0"
P4$Timepoint <- as.numeric(P4$Timepoint)


part_corr_P4 <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  df <- cbind(P4,proteins[,i])
  corr <- pcor(df)
  part_corr_P4[1,i] <- corr[1]$estimate[1,3]
  part_corr_P4[2,i] <- corr[2]$p.value[1,3]
}

rownames(part_corr_P4) <- c("r", "p_val")
colnames(part_corr_P4) <- colnames(proteins)


part_corr_P4 <- t(part_corr_P4)
part_corr_P4 <- data.frame(part_corr_P4)

part_corr_P4$adj_pval <- p.adjust(part_corr_P4$p_val)

part_corr_P4_sig <- part_corr_P4[part_corr_P4$adj_pval < 0.05,]

part_corr_P4_sig <- part_corr_P4_sig[abs(part_corr_P4_sig$r) > 0.5,]

# Estrogen
E2$Timepoint <- md_timepoint$TimePoint
E2$Timepoint[E2$Timepoint=="1st"] <- "1"
E2$Timepoint[E2$Timepoint=="2nd"] <- "2"
E2$Timepoint[E2$Timepoint=="3rd"] <- "3"
E2$Timepoint[E2$Timepoint=="Postpartum"] <- "4"
E2$Timepoint[E2$Timepoint=="Before"] <- "0"
E2$Timepoint[E2$Timepoint=="Non-pregnant"] <- "0"
E2$Timepoint <- as.numeric(E2$Timepoint)

part_corr_E2 <- data.frame(matrix(nrow =2, ncol = length(proteins)), stringsAsFactors = F)

for (i in 1:length(proteins)) {
  df <- cbind(E2,proteins[,i])
  corr <- pcor(df)
  part_corr_E2[1,i] <- corr[1]$estimate[1,3]
  part_corr_E2[2,i] <- corr[2]$p.value[1,3]
}

rownames(part_corr_E2) <- c("r", "p_val")
colnames(part_corr_E2) <- colnames(proteins)


part_corr_E2 <- t(part_corr_E2)
part_corr_E2 <- data.frame(part_corr_E2)

part_corr_E2$adj_pval <- p.adjust(part_corr_E2$p_val)

part_corr_E2_sig <- part_corr_E2[part_corr_E2$adj_pval < 0.05,]

part_corr_E2_sig <- part_corr_E2_sig[abs(part_corr_E2_sig$r) > 0.5,]


# Visualization-----------------------------------------------------------------

# Progesterone
# LIFR
proteins_LIFR <- data.frame(LIF_R=proteins[,"LIF_R", drop=FALSE])

all(rownames(P4)==rownames(proteins_LIFR))

P4_LIFR <- data.frame(P4 = P4$P4,LIFR = proteins_LIFR$LIF_R)
rownames(P4_LIFR) <- rownames(P4)

all(rownames(P4_LIFR)==rownames(md_timepoint))

P4_LIFR$TimePoint <- md_timepoint$TimePoint


ggplot(P4_LIFR, aes(x=LIFR, y=P4, color=TimePoint)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)
  theme_ipsum()

# PDL1
proteins_PDL1 <- data.frame(PDL1=proteins[,"PD_L1", drop=FALSE])
  
all(rownames(P4)==rownames(proteins_PDL1))
  
P4_PDL1 <- data.frame(P4 = P4$P4, PDL1 = proteins_PDL1$PD_L1)
rownames(P4_PDL1) <- rownames(P4)
  
all(rownames(P4_PDL1)==rownames(md_timepoint))
  
P4_PDL1$TimePoint <- md_timepoint$TimePoint
  
  ggplot(P4_PDL1, aes(x=PDL1, y=P4, color=TimePoint)) + 
    geom_point(size=6) +
    geom_smooth(method=lm , color="black", se=FALSE)+
    theme_ipsum()
  
# TWEAK
proteins_TWEAK <- data.frame(TWEAK=proteins[,"TWEAK", drop=FALSE])
  
all(rownames(P4)==rownames(proteins_TWEAK))
  
P4_TWEAK <- data.frame(P4 = P4$P4, TWEAK = proteins_TWEAK$TWEAK)
rownames(P4_TWEAK) <- rownames(P4)
  
all(rownames(P4_TWEAK)==rownames(md_timepoint))
  
P4_TWEAK$TimePoint <- md_timepoint$TimePoint
  
ggplot(P4_TWEAK, aes(x=TWEAK, y=P4, color=TimePoint)) + 
    geom_point(size=6) +
    geom_smooth(method=lm , color="black", se=FALSE)+
    theme_ipsum()
  
# Estrogen
# LIFR
all(rownames(E2)==rownames(proteins_LIFR))

E2_LIFR <- data.frame(E2 = E2$E2,LIFR = proteins_LIFR$LIF_R)
rownames(E2_LIFR) <- rownames(E2)

all(rownames(E2_LIFR)==rownames(md_timepoint))

E2_LIFR$TimePoint <- md_timepoint$TimePoint


ggplot(E2_LIFR, aes(x=LIFR, y=E2, color=TimePoint)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)+
  theme_ipsum()

# PDL1
all(rownames(E2)==rownames(proteins_PDL1))

E2_PDL1 <- data.frame(E2 = E2$E2, PDL1 = proteins_PDL1$PD_L1)
rownames(E2_PDL1) <- rownames(E2)

all(rownames(E2_PDL1)==rownames(md_timepoint))

E2_PDL1$TimePoint <- md_timepoint$TimePoint

ggplot(E2_PDL1, aes(x=PDL1, y=E2, color=TimePoint)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)+
  theme_ipsum()


# TWEAK
all(rownames(E2)==rownames(proteins_TWEAK))

E2_TWEAK <- data.frame(E2 = E2$E2, TWEAK = proteins_TWEAK$TWEAK)
rownames(E2_TWEAK) <- rownames(E2)

all(rownames(E2_TWEAK)==rownames(md_timepoint))

E2_TWEAK$TimePoint <- md_timepoint$TimePoint

ggplot(E2_TWEAK, aes(x=TWEAK, y=E2, color=TimePoint)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE)+
  theme_ipsum()


# Correlation time point by time point------------------------------------------

# P4
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

# E2
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


# Visualize individual proteins per time point

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




# CSF1 3rd trimester E2
proteins_CSF1 <- data.frame(CSF1=proteins[,"CSF_1", drop=FALSE])
CSF1_3rd <- proteins_CSF1[timepoint_3rd,,drop=FALSE]
E2_3rd <- E2[timepoint_3rd,,drop=FALSE]

all(rownames(E2_3rd)==rownames(CSF1_3rd))

E2_CSF1 <- data.frame(E2 = E2_3rd$E2,CSF1 = CSF1_3rd$CSF_1)

ggplot(E2_CSF1, aes(x=CSF1, y=E2,color=md_3rd$Individual)) + 
  geom_point(size=6) +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()

