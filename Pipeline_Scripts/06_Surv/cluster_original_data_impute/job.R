#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to generate and save the cluster labels for the complete data

argv <- as.numeric(commandArgs(T))

######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)

######################################################################
# Initialize folders
homeFolder <- "/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features/cluster_original_data_impute",
                    sep="/")
resultsFolder <- paste(testFolder, "Results", sep="/")
labelsFolder <- paste(resultsFolder, "Labels", sep="/")

######################################################################
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 35
numFeat <- 2000

# Initialize variable parameters
# Data set which will be tested
cancer <- cancerTypes[argv[1]]
# Imputation method
runType <- argv[2]

######################################################################
# Load functions

# Note: some functions depend on variables initialized above!
# lsaImputation:
# -imputedFile, incompleteFile, projectFolder, jvmGBLimit
# iClusterClustering:
# -numCores
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
# Here we take the union of all the cases IDs as we want to use all of our methods on 
# this union 
# Extract all cases which appear in at least one of the data types
unionData <- columnUnion(cases)


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

unionInd <- featureSubsetIndices(unionData)
unionData <- subsetData(unionData, unionInd)
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

unionStat <- rowStatistics(unionData)
unionData <- normalizeData(unionData, unionStat)
######################################################################
### impute on missing data

sampleRows <- FALSE

imputationMethods <- c(knnImputation, llsImputation, lsaImputation,
                       randomImputation)

# Impute the missing data
impute <- function(x) imputationMethods[[runType]](x, sampleRows)
imputedData <- impute(unionData)

clusteringData <- imputedData

#####################################################################
### Select the number of clusters using SNF's cluster size estimate

data <- lapply(clusteringData, t)

# Calculate the distance between samples
distances <- lapply(data, function(x) as.matrix(dist(x)))

# Convert the distances to affinities
affinities <- lapply(distances, affinityMatrix)

# Fuse the affinity matrices
fusedMatrix <- SNF(affinities)

numClusEstimates <- unlist(
  estimateNumberOfClustersGivenGraph(fusedMatrix, 2:10))
numClus <- max(numClusEstimates[c(1,3)])

######################################################################
# Generate and save the cluster labels for the complete data

sampleRows <- FALSE
clusteringMethods <- c(hierarchicalClustering, iClusterClustering,
                       SNFClustering)

for (i in 1:length(clusteringMethods)) {
  cluster <- function(x) clusteringMethods[[i]](x, numClus,
                                                sampleRows)
  labels <- cluster(clusteringData)
  fileName <- paste(paste(c(argv, i), collapse="_"), ".txt", sep="")
  filePath <- paste(labelsFolder, fileName, sep="/")
  write(labels, filePath, ncolumns=length(labels))
}
