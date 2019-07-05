#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to generate and save the cluster labels for the complete data

argv <- as.numeric(commandArgs(T))
######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
######################################################################
# Initialize folders
homeFolder <- "<ENTER PATH TO cluster_complete_data HERE>"
resultsFolder <- paste(homeFolder, "Results", sep="/")
labelsFolder <- paste(resultsFolder, "Labels", sep="/")
######################################################################
# Initialize fixed variables
# Note that cancerTypes and dataTypes should correspond to the name of the data files. 
# For example, if the data type is mrna and the cancer is BRCA, then the data file should be 
# BRCA_mrna.txt
cancerTypes <- "<ENTER NAME OF CANCER HERE>" # for example  c('BRCA', 'KIRC')
clusterSizes <- "<ENTER SIZE OF CLUSTER HERE" # for example c(4,4)
dataTypes <- "<ENTER NAME OF DATA TYPES HERE>" # for example c('methyl', 'mirna', 'mrna')
numCores <- "<ENETER NUMBER OF CORES TO BE USED IN ICLUSTER FUNCTION (RECOMMENDED 35)>" 
numFeat <- "<CHOOSE NUMBER OF FEATURES TO KEEP (RECOMMENDED 2000)>"

# Initialize variable parameters
# Data set which will be tested
cancer <- cancerTypes[argv[1]]
numClus <- clusterSizes[argv[1]]

######################################################################
# Load functions

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

######################################################################
# Load the original data

loadData <- function(dataType, suffix="") {
  fileName <- paste(cancer, "_", dataType, suffix,".txt", sep="")
  filePath <- paste(projectFolder, "Data", fileName, sep="/")
  return(as.matrix(read.delim(filePath)))
}

numViews <- length(dataTypes)
cases <- vector("list", numViews)
controls <- vector("list", numViews)

# Load the biological data
for (v in 1:numViews) {
  cases[[v]] <- loadData(dataTypes[v], "_cases")
  controls[[v]] <- loadData(dataTypes[v], "_controls")
}

######################################################################
# Generate subsets of the cases derived from the set of individuals
# which are present in all data types

# Extract all cases which appear in all of the data types
completeData <- columnIntersection(cases)

######################################################################
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

completeInd <- featureSubsetIndices(completeData)
completeData <- subsetData(completeData, completeInd)

######################################################################
# Normalize the features in the data sets.

rowStatistics <- function(cases) {
  numViews <- length(cases)
  rowStats <- vector("list", numViews)
  
  for (v in 1:numViews) {
    # Calculate the row means and standard deviations
    rowMean <- apply(cases[[v]], 1, mean, na.rm=TRUE)
    rowSd <- apply(cases[[v]], 1, sd, na.rm=TRUE)
    constantInd <- rowSd==0
    rowSd[constantInd] <- 1
    rowStats[[v]] <- list(mean=rowMean, sd=rowSd, ind=constantInd)
  }
  
  return(rowStats)
}

normalizeData <- function(data, stat) {
  for (v in 1:length(data)) {
    data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
    data[[v]] <- data[[v]][!stat[[v]]$ind, ]
  }
  
  return(data)
}

completeStat <- rowStatistics(completeData)
completeData <- normalizeData(completeData, completeStat)

######################################################################
# Generate and save the cluster labels for the complete data

sampleRows <- FALSE
clusteringMethods <- c(hierarchicalClustering, iClusterClustering,
                       SNFClustering)

for (i in 1:length(clusteringMethods)) {
  cluster <- function(x) clusteringMethods[[i]](x, numClus,
                                                sampleRows)
  labels <- cluster(completeData)
  fileName <- paste(paste(c(argv, i), collapse="_"), ".txt", sep="")
  filePath <- paste(labelsFolder, fileName, sep="/")
  write(labels, filePath, ncolumns=length(labels))
}
