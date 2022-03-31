# Principle component analysis--------------------------------------------------

# Author: Sandra Hellberg
# email: sandra.hellberg@liu.se

# Date: 220321

# Script name: PCA


# 1. Importation and formatting the input data
# 2. Filtering and subsetting 
# 3. Calculate delta values for MSvsHP
# 4. Differential expression analysis MS vs HP
# 5. Volcano plots


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


# 2. Filtering and subsetting--------------------------------------------------- 

# subset md file
idx_MS <- grep("MS", md$Individual)
md_MS <- md[idx_MS,]

idx_HP <- grep("HP", md$Individual)
md_HP <- md[idx_HP,]

# subset protein data
idx_MS <- md_MS$Sample_ID
proteins_MS <- proteins[,idx_MS]

idx_HP <- md_HP$Sample_ID
proteins_HP <- proteins[,idx_HP]


# remove HC from the data
idx_HC <- grep("HC", md$Individual)
md_clean <- md[-idx_HC,]

not_HC <- md_clean$Sample_ID
proteins_clean <- proteins[,not_HC]

# 3. Calculate delta values for MSvsHP------------------------------------------

delta <- proteins_clean 
individuals <- unique(md_clean$Individual_ID)
individuals_toremove <- character(0)

timepoint <- unique(md_clean$TimePoint)

for (i in individuals) {
  index_1st <- (md_clean$Individual_ID==i) & (md_clean$TimePoint=="1st")
  if(any(index_1st)) {
    for (t in timepoint){
      index_t <- (md_clean$Individual_ID==i) & (md_clean$TimePoint==t)
      delta[,index_t] <- proteins_clean[,index_t]-proteins_clean[,index_1st]
    }
  } else {individuals_toremove <- c(individuals_toremove,i)}
}

delta <- delta[,md_clean$Individual_ID!=individuals_toremove]
idx_1st <- grep("_1", colnames(delta))
delta <- delta[,-idx_1st]

idx_BP <- grep("_0", colnames(delta))
delta <- delta[,-idx_BP]

idx_tokeep <- colnames(delta)

md_clean <- md_clean[md_clean$Sample_ID %in% colnames(delta),]

# PCA with prcomp---------------------------------------------------------------

# MS
pca_MS <- prcomp(t(proteins_MS)) # to calculate the principle components the methylation matrix must be transposed
autoplot(pca_MS, data=md_MS, colour="TimePoint", size=4) +
  theme_classic() +
  grids(linetype="dashed")

# HP
pca_HP <- prcomp(t(proteins_HP)) # to calculate the principle components the methylation matrix must be transposed
autoplot(pca_HP, data=md_HP, colour="TimePoint", size=4) +
  theme_classic() +
  grids(linetype="dashed")

# delta
pca_delta <- prcomp(t(delta)) # to calculate the principle components the methylation matrix must be transposed
autoplot(pca_delta, data=md_clean, colour="Individual", size=4) +
  theme_classic() +
  grids(linetype="dashed")

# delta only 2nd trimester
idx_2nd <- grep("_2", colnames(delta))
delta_2nd <- delta[,idx_2nd]

idx_tokeep <- colnames(delta_2nd)
md_2nd <- md_clean[md_clean$Sample_ID %in% colnames(delta_2nd),]

pca_2nd <- prcomp(t(delta_2nd)) # to calculate the principle components the methylation matrix must be transposed
autoplot(pca_2nd, data=md_2nd, colour="Individual", size=4) +
  theme_classic() +
  grids(linetype="dashed")

# delta only 3rd trimester
idx_3rd <- grep("_3", colnames(delta))
delta_3rd <- delta[,idx_3rd]

idx_tokeep <- colnames(delta_3rd)
md_3rd <- md_clean[md_clean$Sample_ID %in% colnames(delta_3rd),]

pca_3rd <- prcomp(t(delta_3rd)) # to calculate the principle components the methylation matrix must be transposed
autoplot(pca_3rd, data=md_3rd, colour="Individual", size=4) +
  theme_classic() +
  grids(linetype="dashed")

# delta only PP trimester
idx_PP <- grep("_4", colnames(delta))
delta_PP <- delta[,idx_PP]

idx_tokeep <- colnames(delta_PP)
md_PP <- md_clean[md_clean$Sample_ID %in% colnames(delta_PP),]

pca_PP <- prcomp(t(delta_PP)) # to calculate the principle components the methylation matrix must be transposed
autoplot(pca_PP, data=md_PP, colour="Individual", size=4) +
  theme_classic() +
  grids(linetype="dashed")
