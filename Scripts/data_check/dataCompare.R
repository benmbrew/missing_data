
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
    
    # subset the clinical data so that it corresponds to individuals in the complete data
    completeIds <- colnames(completeData[[1]])
    completeIds <- transformIDFormat(completeIds)
    
    # find the position of the patient IDS in the clinical data 
    clinicalData <- as.data.frame(clinicalData) # not in Daniel's original code 
    clinicalIds <- as.character(clinicalData$bcr_patient_barcode)
    clinicalInd <- match(completeIds, clinicalIds) # returns a vector of positions of (first) matches of its 
    # first argument in its second. Takes length of x with positions of y. NA where x is not in y. So this will
    # be length of complete_ids with position of clinical ids where they match.
    clinicalData <- clinicalData[clinicalInd, ]# now clinical data has ids match with complete data (cases)
    
  } else {
    # Transform patient IDs to the clinical ID format
    transformIDFormat <- function(x) {
      # Keep the first 12 characters
      x <- substr(x, 1, 12)
      # Change each "." to "-"
      x <- gsub(".", "-", x, fixed=TRUE)
      # Make all letters lowercase
      x <- tolower(x)
      
      return(x)
    }
    
    # Extract all cases which appear in at least one of the data types
    unionData <- columnUnion(cases)
    
    # Subset the clinical data so that it corresponds to individuals
    # in the union data
    unionIDs <- colnames(unionData[[1]])
    unionIDs <- transformIDFormat(unionIDs)
    # Find the position of the patient IDs in the clinical data
    clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
    clinicalInd <- match(unionIDs, clinicalIDs)
    clinicalData <- clinicalData[clinicalInd, ]
    
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
  # #   data will still be normalized after imputation (before clustering).
  rowStatistics <- function(cases){
    num_views <- length(cases)
    row_stats <- vector('list', num_views)
    
    for(v in 1:num_views){
      #calculate the row means and std deviations 
      row_mean <- apply(cases[[v]], 1, mean, na.rm = T)
      row_sd <- apply(cases[[v]], 1, sd, na.rm = T)
      constant_ind <- row_sd == 0
      row_sd[constant_ind] <- 1
      row_stats[[v]] <- list(mean = row_mean, sd = row_sd, ind = constant_ind)
    }
    return(row_stats)
  }
  
  normalizeData <- function(data, stat){
    for(v in 1:length(data)) {
      data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
      data[[v]] <- data[[v]][!stat[[v]]$ind, ]
    }
    return(data)
  }
  
  if (complete) {
    completeStat <- rowStatistics(completeData)
    completeData <- normalizeData(completeData, completeStat)
  } else {
    unionStat <- rowStatistics(unionData)
    unionData <- normalizeData(unionData, unionStat)
  }
  
  if (complete) {
    return(list(first = completeData, second = clinicalData))
    
  } else {
    return(list(first = unionData, second = clinicalData))
  }
}

#####################################################################################
# Load complete, union and clinical data.

BRCAComplete <- loadData(cancer = 'BRCA', clinicalDataBRCA, complete = TRUE)
BRCAUnion <- loadData(cancer = 'BRCA', clinicalDataBRCA, complete = FALSE)
clinicalDataCompleteBRCA <- BRCAComplete[[2]] 
clinicalDataUnionBRCA <- BRCAUnion[[2]]
BRCAComplete <- BRCAComplete[[1]]
BRCAUnion <- BRCAUnion[[1]]

KIRCComplete <- loadData(cancer = 'KIRC', clinicalDataKIRC, complete = TRUE)
KIRCUnion <- loadData(cancer = 'KIRC', clinicalDataKIRC, complete = FALSE)
clinicalDataCompleteKIRC <- KIRCComplete[[2]] 
clinicalDataUnionKIRC <- KIRCUnion[[2]]
kircComplete <- KIRCComplete[[1]]
kircUnion <- KIRCUnion[[1]]

LIHCComplete <- loadData(cancer = 'LIHC', clinicalDataLIHC, complete = TRUE)
LIHCUnion <- loadData(cancer = 'LIHC', clinicalDataLIHC, complete = FALSE)
clinicalDataCompleteLIHC <- LIHCComplete[[2]] 
clinicalDataUnionLIHC <- LIHCUnion[[2]]
LIHCComplete <- LIHCComplete[[1]]
LIHCUnion <- LIHCUnion[[1]]

LUADComplete <- loadData(cancer = 'LUAD', clinicalDataLUAD, complete = TRUE)
LUADUnion <- loadData(cancer = 'LUAD', clinicalDataLUAD, complete = FALSE)
clinicalDataCompleteLUAD <- LUADComplete[[2]] 
clinicalDataUnionLUAD <- LUADUnion[[2]]
LUADComplete <- LUADComplete[[1]]
LUADUnion <- LUADUnion[[1]]
