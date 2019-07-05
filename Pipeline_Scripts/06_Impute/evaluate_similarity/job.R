#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to evaluate the how each similarity method affects the
# performance of the clustering methods

argv <- as.numeric(commandArgs(T))

######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(survival)

######################################################################
# Initialize folders
homeFolder <- "/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Impute",
                    "evaluate_similarity", sep="/")
resultsFolder <- paste(testFolder, "Results", sep="/")

######################################################################
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 1
numFeat <- 2000

# Initialize variable parameters
# Data set which will be tested
cancer <- cancerTypes[argv[1]]
cancerInd <- argv[1]

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

# Load the clinical data
clinicalData <- loadData("clin")

######################################################################
# Generate subsets of the cases derived from the set of individuals
# which are present in all data types

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

# Extract all cases which appear in all of the data types
completeData <- columnIntersection(cases)

# Subset the clinical data so that it corresponds to individuals
# in the complete data
completeIDs <- colnames(completeData[[1]])
completeIDs <- transformIDFormat(completeIDs)
# Find the position of the patient IDs in the clinical data
clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
clinicalInd <- match(completeIDs, clinicalIDs)
clinicalData <- clinicalData[clinicalInd, ]

# Generate an incomplete data set by extracting all cases which appear
# in all of the data types and removing some of their samples

set.seed(1)

incompleteIntersection <- generateIncompleteIntersection(cases)
incompleteData <- incompleteIntersection$incomplete
removedData <- incompleteIntersection$removed

# Extract all cases which appear in all of the data types in the
# incomplete data

# Remove samples from the data matrix which only contain NA values.
removeNASamples <- function(data) {
  missingInd <- apply(data, 2, function(x) all(is.na(x)))
  data <- data[, !missingInd]
  return(data)
}
removedNAData <- lapply(incompleteData, removeNASamples)
intersectedData <- columnIntersection(removedNAData)

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

incompleteInd <- featureSubsetIndices(incompleteData)
incompleteData <- subsetData(incompleteData, incompleteInd)
removedData <- subsetData(removedData, incompleteInd)

intersectedInd <- featureSubsetIndices(intersectedData)
intersectedData <- subsetData(intersectedData, intersectedInd)

######################################################################
# Normalize the features in the data sets.
# Normalization is performed before imputation and we expect that the
# data will still be normalized after imputation (before clustering).

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

incompleteStat <- rowStatistics(incompleteData)
incompleteData <- normalizeData(incompleteData, incompleteStat)
removedData <- normalizeData(removedData, incompleteStat)

intersectedStat <- rowStatistics(intersectedData)
intersectedData <- normalizeData(intersectedData, intersectedStat)

######################################################################
# Evaluate the performance of the similarity methods

sampleRows <- FALSE
similarityMethods <- c(selfSimilarity, medianSimilarity,
                       regressionSimilarity)


for (i in 1:length(similarityMethods)) {
  
  similarity <- similarityMethods[[i]]
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    incompleteData<- lapply(incompleteData, t)
  }
  
  # Calculate the distance between samples
  distances <- lapply(incompleteData, function(x) as.matrix(dist(x)))
  
  # Convert the distances to affinities
  affinities <- lapply(distances, affinityMatrix)
  
  # Replace the missing similarities
  affinities <- similarity(affinities)
  
    for(j in 1:length(dataTypes)) {
      write.csv(affinities[[j]], paste0(resultsFolder, '/affinity', '_', argv, '_', i, '_', j, '.csv'))
  }
  
} 
  
