#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to generate and save the cluster labels for the complete data

######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(sva)

######################################################################
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Combat", 
                    "cluster_complete_data", sep="/")
resultsFolder <- paste(testFolder, "Results", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
labelsFolder <- paste(resultsFolder, "Labels", sep="/")

######################################################################
# Initialize fixed variables
cancerTypes <- 'KIRC'
clusterSizes <- 5
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 35
numFeat <- 2000
numViews <- 3


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

#####################################################################
# Generate subsets of the cases derived from the set of individuals
# which are present in all data types

# Extract all cases which appear in all of the data types
completeData <- columnIntersection(cases)

#####################################################################
# subset clinical data by ids in cases
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

###### Run combat function 
for (i in 1:numViews) {
  temp.completeData <- completeData[[i]]
  temp.completeData <- temp.completeData[,completeInd]
  temp.2.completeData <- temp.completeData[rowSums(temp.completeData) != 0,]
  temp.modcombat <- model.matrix(~1, data = clinicalData)
  temp.batch <- clinicalData$gender
  temp_combat = ComBat(dat=temp.2.completeData, batch=temp.batch, mod=temp.modcombat, par.prior=TRUE, prior.plots=FALSE)
  completeData[[i]] <- temp_combat
}

#####################################################################
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
numClus <- 5
sampleRows <- FALSE
clusteringMethods <- c(hierarchicalClustering, iClusterClustering,
                       SNFClustering)

for (i in 1:length(clusteringMethods)) {
  cluster <- function(x) clusteringMethods[[i]](x, numClus,
                                                sampleRows)
  labels <- cluster(completeData)
  fileName <- paste(i, collapse="_", ".txt", sep="")
  filePath <- paste(labelsFolder, fileName, sep="/")
  write(labels, filePath, ncolumns=length(labels))
}
