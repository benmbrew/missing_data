
###############################################################################
# This script will look at the difference between clinical data in the intersection and union.

# Load libraries
library(dplyr)
library(ggplot2)
library(reshape2)

###############################################################################
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Combat",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(projectFolder, 'Scripts/06_Results', sep = '/')

###############################################################################
# Initialize fixed variables
jvmGBLimit <- 8
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 12
numFeat <- 2000

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))
source(paste0(resultsFolder, '/Lib/helpers.R'))

#########################################################################################
# Load the original data
loadData <- function(cancer, clinicalData, complete = FALSE){
  
  loadData <- function(dataType, suffix="") {
    fileName <- paste(cancer, "_", dataType, suffix,".txt", sep="")
    filePath <- paste(projectFolder, "Data", fileName, sep="/")
    return(read.delim(filePath))
  }
  
  numViews <- length(dataTypes)
  cases <- vector("list", numViews)
  controls <- vector("list", numViews)
  
  # Load the biological data
  for (v in 1:numViews) {
    cases[[v]] <- as.matrix(loadData(dataTypes[v], "_cases"))
    controls[[v]] <- as.matrix(loadData(dataTypes[v], "_controls"))
  }
  
  
  
  # transform patient IDs to the clinical ID format 
  transformIDFormat <- function(x){
    x <- substr(x, 1, 12)
    x <- gsub('.', '-', x, fixed = TRUE)
    x <- tolower(x)
    
    return(x)
  }
  
  if (complete) {
    # extract all cases which appear in all of the data types (intersection)
    completeData <- columnIntersection(cases) # now ncol in each matrix is the same and identical to each other. 
    
    
  } else {
  
    
    # Extract all cases which appear in at least one of the data types
    unionData <- columnUnion(cases)
    
    
  }
  
  # Select a subset of features which differ most between cases and
  # controls.
  
  featureSubsetIndices <- function(cases, subsetSize=numFeat) {
    numViews <- length(cases)
    featureSubsetInd <- vector("list", numViews)
    
    for (v in 1:numViews) {
      # Calculate the t-test p-value for each feature, grouped by cases
      # and controls
      numFeatures <- nrow(cases[[v]])
      pval <- sapply(1:numFeatures,
                     function(i) t.test(cases[[v]][i, ],
                                        controls[[v]][i, ])$p.value)
      
      # Subset the data keeping the features with the smallest p-values
      ind <- order(pval)
      featureSubsetInd[[v]] <- ind[1:min(subsetSize, numFeatures)]
    }
    
    return(featureSubsetInd)
  }
  
  subsetData <- function(data, ind) {
    
    for (v in 1:length(data)) {
      data[[v]] <- data[[v]][ind[[v]], ]
    }
    
    return(data)
  }
  
  if (complete){
    completeInd <- featureSubsetIndices(completeData)
    completeData <- subsetData(completeData, completeInd)
  } else {
    unionInd <- featureSubsetIndices(unionData)
    unionData <- subsetData(unionData, unionInd)
  }
  #   
  #   ####################################################################################
  # #   Normalize the features in the data sets.
  # #   Normalization is performed before imputation and we expect that the
#   # #   data will still be normalized after imputation (before clustering).
#   rowStatistics <- function(cases){
#     num_views <- length(cases)
#     row_stats <- vector('list', num_views)
#     
#     for(v in 1:num_views){
#       #calculate the row means and std deviations 
#       row_mean <- apply(cases[[v]], 1, mean, na.rm = T)
#       row_sd <- apply(cases[[v]], 1, sd, na.rm = T)
#       constant_ind <- row_sd == 0
#       row_sd[constant_ind] <- 1
#       row_stats[[v]] <- list(mean = row_mean, sd = row_sd, ind = constant_ind)
#     }
#     return(row_stats)
#   }
#   
#   normalizeData <- function(data, stat){
#     for(v in 1:length(data)) {
#       data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
#       data[[v]] <- data[[v]][!stat[[v]]$ind, ]
#     }
#     return(data)
#   }
#   
#   if (complete) {
#     completeStat <- rowStatistics(completeData)
#     completeData <- normalizeData(completeData, completeStat)
#   } else {
#     unionStat <- rowStatistics(unionData)
#     unionData <- normalizeData(unionData, unionStat)
#   }
  
  if (complete) {
    return(completeData)
    
  } else {
    return(unionData)
  }
}

#####################################################################################
# Load complete, union and clinical data.

BRCAComplete <- loadData(cancer = 'BRCA', complete = TRUE)
BRCAUnion <- loadData(cancer = 'BRCA', complete = FALSE)

KIRCComplete <- loadData(cancer = 'KIRC', complete = TRUE)
KIRCUnion <- loadData(cancer = 'KIRC', complete = FALSE)

LIHCComplete <- loadData(cancer = 'LIHC', complete = TRUE)
LIHCUnion <- loadData(cancer = 'LIHC', complete = FALSE)

LUADComplete <- loadData(cancer = 'LUAD', complete = TRUE)
LUADUnion <- loadData(cancer = 'LUAD', complete = FALSE)

#######################################################################################
# separate into data types 
# BRCA
brca_comp_methyl <- BRCAComplete[[1]]
brca_comp_mirna <- BRCAComplete[[2]]
brca_comp_mrna <- BRCAComplete[[3]]

brca_union_methyl <- BRCAUnion[[1]]
brca_union_mirna <- BRCAUnion[[2]]
brca_union_mrna <- BRCAUnion[[3]]

# kirc
kirc_comp_methyl <- KIRCComplete[[1]]
kirc_comp_mirna <- KIRCComplete[[2]]
kirc_comp_mrna <- KIRCComplete[[3]]

kirc_union_methyl <- KIRCUnion[[1]]
kirc_union_mirna <- KIRCUnion[[2]]
kirc_union_mrna <- KIRCUnion[[3]]

# lihc
lihc_comp_methyl <- LIHCComplete[[1]]
lihc_comp_mirna <- LIHCComplete[[2]]
lihc_comp_mrna <- LIHCComplete[[3]]

lihc_union_methyl <- LIHCUnion[[1]]
lihc_union_mirna <- LIHCUnion[[2]]
lihc_union_mrna <- LIHCUnion[[3]]

# luad
luad_comp_methyl <- LUADComplete[[1]]
luad_comp_mirna <- LUADComplete[[2]]
luad_comp_mrna <- LUADComplete[[3]]

luad_union_methyl <- LUADUnion[[1]]
luad_union_mirna <- LUADUnion[[2]]
luad_union_mrna <- LUADUnion[[3]]

#######################################################################################
# extract samples in union that are not in intersection for each data type 

# BRCA
brca_not_methyl <- brca_union_methyl[,!colnames(brca_union_methyl) %in% colnames(brca_comp_methyl)]
brca_not_mirna <- brca_union_mirna[,!colnames(brca_union_mirna) %in% colnames(brca_comp_mirna)]
brca_not_mrna <- brca_union_mrna[,!colnames(brca_union_mrna) %in% colnames(brca_comp_mrna)]

# KIRC
kirc_not_methyl <- kirc_union_methyl[,!colnames(kirc_union_methyl) %in% colnames(kirc_comp_methyl)]
kirc_not_mirna <- kirc_union_mirna[,!colnames(kirc_union_mirna) %in% colnames(kirc_comp_mirna)]
kirc_not_mrna <- kirc_union_mrna[,!colnames(kirc_union_mrna) %in% colnames(kirc_comp_mrna)]

# lihc
lihc_not_methyl <- lihc_union_methyl[,!colnames(lihc_union_methyl) %in% colnames(lihc_comp_methyl)]
lihc_not_mirna <- lihc_union_mirna[,!colnames(lihc_union_mirna) %in% colnames(lihc_comp_mirna)]
lihc_not_mrna <- lihc_union_mrna[,!colnames(lihc_union_mrna) %in% colnames(lihc_comp_mrna)]

# luad
luad_not_methyl <- luad_union_methyl[,!colnames(luad_union_methyl) %in% colnames(luad_comp_methyl)]
luad_not_mirna <- luad_union_mirna[,!colnames(luad_union_mirna) %in% colnames(luad_comp_mirna)]
luad_not_mrna <- luad_union_mrna[,!colnames(luad_union_mrna) %in% colnames(luad_comp_mrna)]


###############################################################################################################
################# Compare brca
par(mfrow = c(1,2))
# methyl
brca_comp_methyl_mean <- apply(brca_comp_methyl, 2, mean)
hist(brca_comp_methyl_mean)
brca_not_methyl_mean <- apply(brca_not_methyl, 2, mean)
hist(brca_not_methyl_mean)

# mirna
brca_comp_mirna_mean <- apply(brca_comp_mirna, 2, mean)
hist(brca_comp_mirna_mean)
brca_not_mirna_mean <- apply(brca_not_mirna, 2, mean)
hist(brca_not_mirna_mean)

# mrna
brca_comp_mrna_mean <- apply(brca_comp_mrna, 2, mean)
hist(brca_comp_mrna_mean)
brca_not_mrna_mean <- apply(brca_not_mrna, 2, mean)
hist(brca_not_mrna_mean)

############### Compare kirc

# methyl
kirc_comp_methyl_mean <- apply(kirc_comp_methyl, 2, mean)
hist(kirc_comp_methyl_mean)
kirc_not_methyl_mean <- apply(kirc_not_methyl, 2, mean)
hist(kirc_not_methyl_mean)

# mirna
kirc_comp_mirna_mean <- apply(kirc_comp_mirna, 2, mean)
hist(kirc_comp_mirna_mean)
kirc_not_mirna_mean <- apply(kirc_not_mirna, 2, mean)
hist(kirc_not_mirna_mean)

# mrna
kirc_comp_mrna_mean <- apply(kirc_comp_mrna, 2, mean)
hist(kirc_comp_mrna_mean)
kirc_not_mrna_mean <- apply(kirc_not_mrna, 2, mean)
hist(kirc_not_mrna_mean)

############### Compare lihc

# methyl
lihc_comp_methyl_mean <- apply(lihc_comp_methyl, 2, mean)
hist(lihc_comp_methyl_mean)
lihc_not_methyl_mean <- apply(lihc_not_methyl, 2, mean)
hist(lihc_not_methyl_mean)

# mirna
lihc_comp_mirna_mean <- apply(lihc_comp_mirna, 2, mean)
hist(lihc_comp_mirna_mean)
lihc_not_mirna_mean <- apply(lihc_not_mirna, 2, mean)
hist(lihc_not_mirna_mean)

# mrna
lihc_comp_mrna_mean <- apply(lihc_comp_mrna, 2, mean)
hist(lihc_comp_mrna_mean)
lihc_not_mrna_mean <- apply(lihc_not_mrna, 2, mean)
hist(lihc_not_mrna_mean)

############### Compare luad

# methyl
luad_comp_methyl_mean <- apply(luad_comp_methyl, 2, mean)
hist(luad_comp_methyl_mean)
luad_not_methyl_mean <- apply(luad_not_methyl, 2, mean)
hist(luad_not_methyl_mean)

# mirna
luad_comp_mirna_mean <- apply(luad_comp_mirna, 2, mean)
hist(luad_comp_mirna_mean)
luad_not_mirna_mean <- apply(luad_not_mirna, 2, mean)
hist(luad_not_mirna_mean)

# mrna
luad_comp_mrna_mean <- apply(luad_comp_mrna, 2, mean)
hist(luad_comp_mrna_mean)
luad_not_mrna_mean <- apply(luad_not_mrna, 2, mean)
hist(luad_not_mrna_mean)

