#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to evaluate the how each similarity method affects the
# performance of the clustering methods

argv <- as.numeric(commandArgs(T))

######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(survival)
library(sva)

######################################################################
# Initialize folders
homeFolder <- "/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "missing_data_combat",
                    "evaluate_similarity", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(testFolder, "Results", sep="/")

######################################################################
# Initialize fixed variables
cancer <- 'KIRC'
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 1
numFeat <- 2000
clusters <- c(1,2,3,4)

# Seed for generating incomplete data
seed <- argv[1]

# Store the output in subfolders
resultsFile <- paste(paste(argv, collapse="_"), ".txt", sep="")
constructPath <- function(intermediateFolders, parent=resultsFolder,
                          file=resultsFile) {
  paste(parent, intermediateFolders, file, sep="/")
}
similarityFile <- constructPath("Similarity")

#####################################################################
# Load in raw clinical data for kirc
transposeDataFrame <- function(df, colnamesInd=1) {
  variableNames <- df[, colnamesInd]
  df <- as.data.frame(t(df[, -colnamesInd]))
  colnames(df) <- variableNames
  df
}


extractRelevantColumns <- function(data) {
  # List of features used for survival analysis
  features <- c("admin.batch_number",
                "patient.bcr_patient_barcode",
                "patient.bcr_patient_uuid",
                "patient.days_to_death",
                "patient.days_to_last_followup",
                "patient.days_to_last_known_alive",
                "patient.vital_status",
                "patient.gender")
  patientFeatures <- paste("patient", features, sep=".")
  
  # Add missing features to the data
  missingFeaturesInd <- !(features %in% colnames(data))
  data[features[missingFeaturesInd]] <- NA
  
  # Extract and rename the relevant columns
  data <- data[features]
  colnames(data) <- features
  
  return(data)
}


loadClinData <- function(cancer = 'KIRC') {
  #   processingResult <- dataFolder
  fileSuffix <- "clinical.txt"
  
  # Load the data
  data <- NULL
  fileName <- paste(cancer, fileSuffix, sep="_")
  filePath <- paste(dataFolder,fileName,
                    sep="/")
  if (file.exists(filePath)) {
    data <- read.delim(filePath)
    
    # Transpose the data
    data <- transposeDataFrame(data)
  }
  
  return(data)
}

# Load clinical data with batch number
cancer <- 'KIRC'
# Process the clinical data
clinicalData <- loadClinData(cancer)
clinicalData <- extractRelevantColumns(clinicalData)

## remove patient from column names
features <- colnames(clinicalData)
split <- strsplit(features, '.', fixed = TRUE)
keepSplit <- lapply(split, function(x) x[length(x)])
features <- unlist(keepSplit)
colnames(clinicalData) <- features

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
####################################################################

# Extract all cases which appear in all of the data types
completeData <- columnIntersection(cases)

######################################################################
# Generate subsets of the cases derived from the set of individuals
# which are present in all data types

transformIDFormat <- function(x) {
  # Keep the first 12 characters
  x <- substr(x, 1, 12)
  # Change each "." to "-"
  x <- gsub(".", "-", x, fixed=TRUE)
  # Make all letters lowercase
  x <- tolower(x)
  
  return(x)
}


completeIDs <- colnames(completeData[[1]])
completeIDs <- transformIDFormat(completeIDs)
# Find the position of the patient IDs in the clinical data
clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
clinicalInd <- match(completeIDs, clinicalIDs)
clinicalData <- clinicalData[clinicalInd, ]
# Subset clinical data and completeData by clinicalData
clinicalData <- clinicalData[rowSums(is.na(clinicalData)) < ncol(clinicalData),]
clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
completeInd <- match(clinicalIDs, completeIDs)

clinicalData$days_to_death <- as.numeric(clinicalData$days_to_death)
clinicalData$days_to_last_followup <- as.numeric(clinicalData$days_to_last_followup)
###### subset data 
for (i in 1:numViews) {
  temp.completeData <- completeData[[i]]
  temp.completeData <- temp.completeData[,completeInd]
  temp.2.completeData <- temp.completeData[rowSums(temp.completeData) != 0,]
  completeData[[i]] <- temp.2.completeData
}

# run combat on methylation
temp.dat <- completeData[[1]]
temp.modcombat <- model.matrix(~1, data = clinicalData)
temp.batch <- clinicalData$gender
temp_combat = ComBat(dat=temp.dat, batch=temp.batch, mod=temp.modcombat, par.prior=TRUE, prior.plots=FALSE)
completeData[[1]] <- temp_combat


# Generate an incomplete data set by extracting all cases which appear
# in all of the data types and removing some of their samples

set.seed(seed)

incompleteIntersection <- generateIncompleteIntersectionCombat(cases)
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

# Read in the SNF cluster labels for the complete data
for (i in 1:length(clusters)) { 
  fileName <- paste(i, ".txt", sep="")
  filePath <- paste(testFolder, "../cluster_complete_data",
                    "Results/Labels", fileName, sep="/")
  completeLabels <- scan(filePath)
  
  # Save the results of clustering the incomplete and intersected data
  for (j in 1:length(similarityMethods)) {
    similarity <- similarityMethods[[j]]
    
    data <- incompleteData
    tcgaID <- function(x) colnames(x[[1]])
    dataInd <- match(tcgaID(data), tcgaID(completeData))
    similarityResults <- evaluateSimilarity(similarity, data,
                                            completeLabels[dataInd],
                                            clinicalData[dataInd, ],
                                            sampleRows)
    writeResults(c(argv, j, i, similarityResults), similarityFile)
  }
}
