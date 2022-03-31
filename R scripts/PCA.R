# Principle component analysis--------------------------------------------------

# Author: Sandra Hellberg
# email: sandra.hellberg@liu.se

# Date: 220321

# 1. Importation and formatting the input data
# 2. Filtering and subsetting 
# 3. Calculate delta values for MSvsHP
# 4. Differential expression analysis MS vs HP
# 5. Volcano plots


# Load packages-----------------------------------------------------------------
pack_R <- c("stats", "ggfortify", "ggplot2", "dplyr", "ggpubr")

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

# 2. Filtering and subsetting--------------------------------------------------- 

# remove healthy controls (HC) from the data prior to downstream analysis
idx_HC <- grep("HC", md$Individual)
md_clean <- md[-idx_HC,]

not_HC <- md_clean$Sample_ID
proteins_clean <- proteins[,not_HC]

all(md_clean$Sample_ID == colnames(proteins_clean))

# subset MS and healthy pregnant (HP)
idx_MS <- grep("MS", md_clean$Individual)
md_MS <- md_clean[idx_MS,]

idx_HP <- grep("HP", md_clean$Individual)
md_HP <- md_clean[idx_HP,]

# subset protein data
idx_MS <- md_MS$Sample_ID
proteins_MS <- proteins_clean[,idx_MS]

idx_HP <- md_HP$Sample_ID
proteins_HP <- proteins_clean[,idx_HP]

# 3. Calculate delta values for MSvsHP------------------------------------------
# The comparison between MS and HP was done using the delta values where the
# 1st trimester time point was substracted from the other time points prior to 
# analysis, i.e 2nd-1st, 3rd-1st, PP-1st. 

delta <- proteins_clean # create a new data frame of the protein values
individuals <- unique(md_clean$Individual_ID) # identify all unique individuals present in the data
individuals_toremove <- character(0) # initialize vector for individuals without 1st time point
timepoint <- unique(md_clean$TimePoint)

# For loop for creating delta values for each individual per time point
for (i in individuals) {
  index_1st <- (md_clean$Individual_ID==i) & (md_clean$TimePoint=="1st")
  if(any(index_1st)) {
    for (t in timepoint){
      index_t <- (md_clean$Individual_ID==i) & (md_clean$TimePoint==t)
      delta[,index_t] <- proteins_clean[,index_t]-proteins_clean[,index_1st] # substract 1st trimester from 2nd, 3rd and PP
    }
  } else {individuals_toremove <- c(individuals_toremove,i)} 
}

delta <- delta[,md_clean$Individual_ID!=individuals_toremove] # remove individual without first trimester time point
idx_1st <- grep("_1", colnames(delta))
delta <- delta[,-idx_1st] # remove the 1st trimester time point

idx_BP <- grep("_0", colnames(delta)) # remove the before pregnancy samples
delta <- delta[,-idx_BP]

md_clean <- md_clean[md_clean$Sample_ID %in% colnames(delta),] #subset the meta data to only include samples with delta values

all(md_clean$Sample_ID==colnames(delta))

# PCA with prcomp---------------------------------------------------------------

# Delta values with all time points, colored by time point
pca_delta <- prcomp(t(delta)) # transpose protein matrix
autoplot(pca_delta, data=md_clean, colour="TimePoint", size=10) +
  theme_classic() +
  grids(linetype="dashed")

# delta only 2nd trimester
idx_2nd <- grep("_2", colnames(delta))
delta_2nd <- delta[,idx_2nd]

idx_tokeep <- colnames(delta_2nd)
md_2nd <- md_clean[md_clean$Sample_ID %in% colnames(delta_2nd),]

all(md_2nd$Sample_ID==colnames(delta_2nd))

pca_2nd <- prcomp(t(delta_2nd)) 
autoplot(pca_2nd, data=md_2nd, colour="Individual", size=10) +
  theme_classic() +
  grids(linetype="dashed")

# delta only 3rd trimester
idx_3rd <- grep("_3", colnames(delta))
delta_3rd <- delta[,idx_3rd]

idx_tokeep <- colnames(delta_3rd)
md_3rd <- md_clean[md_clean$Sample_ID %in% colnames(delta_3rd),]

pca_3rd <- prcomp(t(delta_3rd)) 
autoplot(pca_3rd, data=md_3rd, colour="Individual", size=10) +
  theme_classic() +
  grids(linetype="dashed")

all(md_3rd$Sample_ID==colnames(delta_3rd))


# delta only PP trimester
idx_PP <- grep("_4", colnames(delta))
delta_PP <- delta[,idx_PP]

idx_tokeep <- colnames(delta_PP)
md_PP <- md_clean[md_clean$Sample_ID %in% colnames(delta_PP),]

pca_PP <- prcomp(t(delta_PP)) 
autoplot(pca_PP, data=md_PP, colour="Individual", size=10) +
  theme_classic() +
  grids(linetype="dashed")

all(md_PP$Sample_ID==colnames(delta_PP))

# End of script-----------------------------------------------------------------