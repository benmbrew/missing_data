#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to evaluate the how each imputation method affects the
# performance of the clustering methods on the original data

argv <- as.numeric(commandArgs(T))

######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)

######################################################################
# Initialize folders
homeFolder <- "/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features_Dup",
                    "evaluate_original_imputation", sep="/")
resultsFolder <- paste(testFolder, "Results", sep="/")

######################################################################
# Initialize fixed variables
jvmGBLimit <- 8
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 12
numFeat <- 2000

# Initialize variable parameters
# Data set which will be tested
cancer <- cancerTypes[argv[1]]
cancerInd <- argv[1]
# Imputation method
runType <- argv[2]

# Store the output in subfolders
resultsFile <- paste(paste(argv, collapse="_"), ".txt", sep="")
constructPath <- function(intermediateFolders, parent=resultsFolder,
                          file=resultsFile) {
  paste(parent, intermediateFolders, file, sep="/")
}
clusteringFile <- constructPath("Clustering")

incompleteFile <- constructPath("Incomplete")
imputedFile <- constructPath("Imputed")

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
# Generate the union and intersected data from the original data

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

transformUnionId <- function(x) {
  # Keep the first 12 characters
  x <- substr(x, 1, 12)
  
  return(x)
}

unionData <- columnUnion(cases)

# Subset the clinical data so that it corresponds to individuals
# in the union data

for (i in 1:3) {
  temp.data  <- unionData[[i]]
  temp.names <- transformUnionId(colnames(temp.data))
  colnames(temp.data) <- temp.names
  temp.data <- temp.data[, !duplicated(colnames(temp.data))]
  unionData[[i]] <- temp.data
}

# Subset the clinical data so that it corresponds to individuals
# in the union data
unionIDs <- colnames(unionData[[1]])
unionIDs <- transformIDFormat(unionIDs)
# Find the position of the patient IDs in the clinical data
clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
clinicalInd <- match(unionIDs, clinicalIDs)
clinicalData <- clinicalData[clinicalInd, ]

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

unionStat <- rowStatistics(unionData)
unionData <- normalizeData(unionData, unionStat)

######################################################################
# Evaluate the performance of the methods

sampleRows <- FALSE
clusteringMethods <- c(hierarchicalClustering, iClusterClustering,
                       SNFClustering)
imputationMethods <- c(knnImputation, llsImputation, lsaImputation,
                       randomImputation)

# Impute the missing data
impute <- function(x) imputationMethods[[runType]](x, sampleRows)
imputedData <- impute(unionData)

clusteringData <- list(imputedData)

# Save the results of clustering the imputed and intersected data
for (i in 1:length(clusteringMethods)) {
  # Read in the cluster labels for the complete data
  fileName <- paste(cancerInd, "_", i, ".txt", sep="")
  filePath <- paste(testFolder, "../cluster_complete_data",
                    "Results/Labels", fileName, sep="/")
  completeLabels <- scan(filePath)
  
  # Extend the completeLabels
  nSamples <- ncol(unionData[[1]])
  repTimes <- ceiling(nSamples/length(completeLabels))
  completeLabels <- rep.int(completeLabels, repTimes)[1:nSamples]
  
  # Initialize clustering function
  numClus <- length(unique(completeLabels))
  cluster <- function(x) clusteringMethods[[i]](x, numClus,
                                                sampleRows)

  for (j in 1:length(clusteringData)) {
    data <- clusteringData[[j]]
    tcgaID <- function(x) colnames(x[[1]])
    dataInd <- match(tcgaID(data), tcgaID(imputedData))
    
    clusteringResults <- evaluateClustering(cluster, data,
                                            completeLabels[dataInd],
                                            clinicalData[dataInd, ])
    writeResults(c(argv, i, j, clusteringResults), clusteringFile)
  }
}
