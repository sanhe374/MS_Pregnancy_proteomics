# Differential expression analysis with limma-----------------------------------

# Author: Georgia Papapavlou and Sandra Hellberg
# email: georgia.papapavlou.linghed@liu.se ; sandra.hellberg@liu.se

# Date: 211207

# Script description: This is a script to perform differential expression analysis to identify
# differentially expressed proteins between trimesters during pregnancy and patients and controls. Results are
# illustrated using Volcano plots. 


# 1. Importation and formatting the input data
# 2. Filtering and subsetting 
# 3. Differential expression analysis across pregnancy
# 4. Differential expression analysis MS vs HP
# 5. Volcano plots

pack_R <- c("BiocParallel", "sva", "tidyverse", "readxl", "assertthat", "ggplot2", "ggfortify", "patchwork",
            "limma", "edgeR", "RColorBrewer", "statmod", "Glimma", "dplyr", "stringr", "tidyr", "ggrepel",
            "EnchancedVolcano")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}


# 1. Import and format input data-----------------------------------------------
md <- read_excel("Md_Olink_newID.xlsx") # read in meta data file
md <- data.frame(md) # convert to a data frame

proteins <- read_excel("proteins_Olink_newID.xlsx") # read in protein data
proteins <- data.frame(proteins)

rownames(proteins) <- proteins$mediator_name # set protein name to row names
proteins <- proteins[,-1]    
proteins_numeric <- sapply(proteins, as.numeric) # convert to numeric
rownames(proteins_numeric) <- rownames(proteins)
proteins <- as.matrix(proteins_numeric)

stopifnot(colnames(proteins) == md$Sample_ID)

# Remove P4 and E2 from the protein data
row_names_df_to_remove<-c("P4","E2")
proteins <- proteins[!(row.names(proteins) %in% row_names_df_to_remove),]

# 2. Filtering and subsetting--------------------------------------------------- 

# Exclude the healthy non-pregnant controls for downstream analysis
is_HC <- which(is.element(md$Individual, "HC"))
proteins <- proteins[, -is_HC]
md <- md[-is_HC, ] # remove the healthy controls

# Subset the data into MS and healthy pregnant
is_MS <- is.element(md$Status, "MS")
proteins_MS <- proteins[, is_MS]
md_MS <- md[is_MS, ]

stopifnot(ncol(proteins_MS) == nrow(md_MS))

is_Healthy <- is.element(md$Status, "Healthy")
proteins_Healthy <- proteins[, is_Healthy]
md_Healthy <- md[is_Healthy, ]

stopifnot(ncol(proteins_Healthy) == nrow(md_Healthy))

# 3. Differential expression analysis across pregnancy--------------------------

# MS
design <- model.matrix(~ 0 + TimePoint, md_MS)
corfit <- duplicateCorrelation(proteins_MS, design, block=md_MS$Individual_ID) # duplicateCorrelation to block the individual effect

contr_matrix <- makeContrasts(
  TimePoint_1stvsBP = TimePoint1st-TimePointBefore,
  TimePoint_2ndvs1st = TimePoint2nd-TimePoint1st,
  TimePoint_3rdvs1st = TimePoint3rd-TimePoint1st,
  TimePoint_PPvs1st = TimePointPostpartum-TimePoint1st,
  TimePoint_2ndvsBP = TimePoint2nd-TimePointBefore,
  TimePoint_3rdvs2nd = TimePoint3rd-TimePoint2nd,
  TimePoint_PPvs2nd = TimePointPostpartum-TimePoint2nd,
  TimePoint_3rdvsBP = TimePoint3rd-TimePointBefore,
  TimePoint_PPvs3rd = TimePointPostpartum-TimePoint3rd,
  TimePoint_PPvsBP = TimePointPostpartum-TimePointBefore,
  levels=colnames(design)
)

fit <- lmFit(proteins_MS , design, block=md_MS$Individual_ID, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
plotSA(fit)
summary(decideTests(fit))

# Extract data for ALL proteins at each comparison
tab_1stvsBP_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_1stvsBP",number=Inf) 
tab_3rdvs1st_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_3rdvs1st",number=Inf)
tab_PPvs1st_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_PPvs1st",number=Inf)
tab_2ndvsBP_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_2ndvsBP",number=Inf)
tab_3rdvs2nd_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_3rdvs2nd",number=Inf)
tab_PPvs2nd_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_PPvs2nd",number=Inf)
tab_3rdvsBP_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_3rdvsBP",number=Inf)
tab_PPvs3rd_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_PPvs3rd",number=Inf)
tab_PPvsBP_MS <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_PPvsBP",number=Inf)

# Extract significant proteins and write to table
names <- c("tab_1stvsBP_MS","tab_3rdvs1st_MS","tab_PPvs1st_MS", "tab_2ndvsBP_MS",
           "tab_3rdvs2nd_MS", "tab_PPvs2nd_MS", "tab_3rdvsBP_MS", "tab_PPvs3rd_MS",
           "tab_PPvsBP_MS")
output <- names %>%
  str_remove("tab_")%>%
  paste0(".txt")

input <- list(tab_1stvsBP_MS, tab_3rdvs1st_MS, tab_PPvs1st_MS, tab_2ndvsBP_MS,
              tab_3rdvs2nd_MS, tab_PPvs2nd_MS, tab_3rdvsBP_MS, tab_PPvs3rd_MS,
              tab_PPvsBP_MS)

for (i in 1:length(input)) {
  sig <- input[[i]][input[[i]]$adj.P.Val < 0.05,]
  write.table(sig,output[i]) 
}

# Healthy pregnant
design <- model.matrix(~ 0 + TimePoint, md_Healthy)
corfit <- duplicateCorrelation(proteins_Healthy, design, block=md_Healthy$Individual_ID)

contr_matrix <- makeContrasts(
  TimePoint_2ndvs1st = TimePoint2nd-TimePoint1st,
  TimePoint_3rdvs1st = TimePoint3rd-TimePoint1st,
  TimePoint_PPvs1st = TimePointPostpartum-TimePoint1st,
  TimePoint_3rdvs2nd = TimePoint3rd-TimePoint2nd,
  TimePoint_PPvs2nd = TimePointPostpartum-TimePoint2nd,
  TimePoint_PPvs3rd = TimePointPostpartum-TimePoint3rd,
  levels=colnames(design)
)


fit <- lmFit(proteins_Healthy , design, block=md_Healthy$Individual_ID, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
plotSA(fit)
summary(decideTests(fit))

# Extract ALL proteins at each time-point
tab_2ndvs1st_HP <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_2ndvs1st",number=Inf)
tab_3rdvs1st_HP  <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_3rdvs1st",number=Inf)
tab_PPvs1st_HP  <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_PPvs1st",number=Inf)
tab_3rdvs2nd_HP  <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_3rdvs2nd",number=Inf)
tab_PPvs2nd_HP  <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_PPvs2nd",number=Inf)
tab_PPvs3rd_HP  <- topTable(fit,adjust.method="BH", p.value = Inf, coef="TimePoint_PPvs3rd",number=Inf)

# Extract significant proteins and write to table
names <- c("tab_2ndvs1st_HP","tab_3rdvs1st_HP","tab_PPvs1st_HP", "tab_3rdvs2nd_HP",
           "tab_PPvs2nd_HP", "tab_PPvs3rd_HP") 

output <- names %>%
  str_remove("tab_")%>%
  paste0(".txt")

input <- list(tab_2ndvs1st_HP, tab_3rdvs1st_HP, tab_PPvs1st_HP, tab_3rdvs2nd_HP,
              tab_PPvs2nd_HP, tab_PPvs3rd_HP)

for (i in 1:length(input)) {
  sig <- input[[i]][input[[i]]$adj.P.Val < 0.05,]
  write.table(sig,output[i]) 
}



# 4. Differential expression analysis MSvsHP------------------------------------
# All comparisons are made after correcting for the first trimester timepoint 

design <- model.matrix(~ 0 + Sample_group, md)
corfit <- duplicateCorrelation(proteins, design, block=md$Individual_ID) #block the individual effect
fit <- lmFit(proteins, design, block=md$Individual_ID, correlation=corfit$consensus)

contr_matrix <- makeContrasts(
  MS2ndvsHP2nd = (Sample_groupMS_2nd - Sample_groupMS_1st) - (Sample_groupHP_2nd - Sample_groupHP_1st),
  MS3rdvsHP3rd = (Sample_groupMS_3rd - Sample_groupMS_1st) - (Sample_groupHP_3rd - Sample_groupHP_1st),
  MSPPvsHPPP = (Sample_groupMS_PP - Sample_groupMS_1st) - (Sample_groupHP_PP - Sample_groupHP_1st),
  levels=colnames(design)
)

fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

# Extract nominally significant proteins for each comparison
DEG_2nd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS2ndvsHP2nd")
DEG_3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS3rdvsHP3rd")
DEG_PP <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MSPPvsHPPP")



# 5. Volcano plots--------------------------------------------------------------

# MS 
# 3rdvst1st
data_in <- tab_3rdvs1st_MS # MS 3rdvs1st
data_in$diffexpressed <- "NO" # create an empty column in the data frame
data_in$diffexpressed[data_in$logFC > 0 & data_in$adj.P.Val < 0.05] <- "UP"
data_in$diffexpressed[data_in$logFC < 0 & data_in$adj.P.Val < 0.05] <- "DOWN"

p <- data_in %>%
  ggplot(mapping=aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed))+
  geom_point(size = 2)+
  theme_minimal()

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="grey", linetype =2)+
  geom_vline(xintercept= c(0), col="grey", linetype =2)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

p3 <- p2 +
  scale_color_manual(values=mycolors)

p3

# MS 
# PPvs3rd

data_in <- tab_PPvs3rd_MS
data_in$diffexpressed <- "NO" # create an empty column in the data frame
data_in$diffexpressed[data_in$logFC > 0 & data_in$adj.P.Val < 0.05] <- "UP"
data_in$diffexpressed[data_in$logFC < 0 & data_in$adj.P.Val < 0.05] <- "DOWN"

p <- data_in %>%
  ggplot(mapping=aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed))+
  geom_point(size = 2)+
  theme_minimal()

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="grey", linetype =2)+
  ##geom_text_repel(data=head(data_in, 20), aes(label=protein))
  geom_vline(xintercept= c(0), col="grey", linetype =2)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

p3 <- p2+
  scale_color_manual(values=mycolors)

p3

# HP 
# 3rdvst1st
data_in <- tab_3rdvs1st_HP
data_in$diffexpressed <- "NO" # create an empty column in the data frame
data_in$diffexpressed[data_in$logFC > 0 & data_in$adj.P.Val < 0.05] <- "UP"
data_in$diffexpressed[data_in$logFC < 0 & data_in$adj.P.Val < 0.05] <- "DOWN"

p <- data_in %>%
  ggplot(mapping=aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed))+
  geom_point(size = 2)+
  theme_minimal()

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="grey", linetype =2)+
  geom_vline(xintercept= c(0), col="grey", linetype =2)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

p3 <- p2+
  scale_color_manual(values=mycolors)

p3


# HP 
# PPvs3rd
data_in <- tab_PPvs3rd_HP
data_in$diffexpressed <- "NO" # create an empty column in the data frame
data_in$diffexpressed[data_in$logFC > 0 & data_in$adj.P.Val < 0.05] <- "UP"
data_in$diffexpressed[data_in$logFC < 0 & data_in$adj.P.Val < 0.05] <- "DOWN"

p <- data_in %>%
  ggplot(mapping=aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed))+
  geom_point(size = 2)+
  theme_minimal()

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="grey", linetype =2)+
  ##geom_text_repel(data=head(data_in, 20), aes(label=protein))
  geom_vline(xintercept= c(0), col="grey", linetype =2)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

p3 <- p2+
  scale_color_manual(values=mycolors)

p3


# End of script-----------------------------------------------------------------

