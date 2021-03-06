#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to evaluate the how each similarity method affects the
# performance of the clustering methods on the original data

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
                    "06_Int",
                    "evaluate_original_similarity", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
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

# Store the output in subfolders
resultsFile <- paste(paste(argv, collapse="_"), ".txt", sep="")
constructPath <- function(intermediateFolders, parent=resultsFolder,
                          file=resultsFile) {
  paste(parent, intermediateFolders, file, sep="/")
}
similarityFile <- constructPath("Similarity")

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

## Use combat if cancer is KIRC 
if (cancer == 'KIRC') {
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
  
  # Process the clinical data
  clinicalData <- loadClinData(cancer)
  clinicalData <- extractRelevantColumns(clinicalData)
  
  ## remove patient from column names
  features <- colnames(clinicalData)
  split <- strsplit(features, '.', fixed = TRUE)
  keepSplit <- lapply(split, function(x) x[length(x)])
  features <- unlist(keepSplit)
  colnames(clinicalData) <- features
  
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
  
  
  # Subset clinical data and completeData by clinicalData
  clinicalData <- clinicalData[rowSums(is.na(clinicalData)) < ncol(clinicalData),]
  clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
  unionInd <- match(clinicalIDs, unionIDs)
  
  clinicalData$days_to_death <- as.numeric(clinicalData$days_to_death)
  clinicalData$days_to_last_followup <- as.numeric(clinicalData$days_to_last_followup)
  ###### Run combat function 
  for (i in 1:numViews) {
    temp.unionData <- unionData[[i]]
    temp.unionData <- temp.unionData[,unionInd]
    temp.2.unionData <- temp.unionData[rowSums(temp.unionData, na.rm = T) != 0,]
    temp.modcombat <- model.matrix(~1, data = clinicalData)
    temp.batch <- clinicalData$gender
    temp_combat = ComBat(dat=temp.2.unionData, batch=temp.batch, mod=temp.modcombat, par.prior=TRUE, prior.plots=FALSE)
    unionData[[i]] <- temp_combat
  }
  
  
  # Subset the clinical data so that it corresponds to individuals
  # in the union data
  
  # Extract all cases which appear in all of the data types
  intersectedData <- columnIntersection(cases)
  
  for (i in 1:3) {
    temp.data  <- intersectedData[[i]]
    temp.names <- transformUnionId(colnames(temp.data))
    colnames(temp.data) <- temp.names
    temp.data <- temp.data[, !duplicated(colnames(temp.data))]
    intersectedData[[i]] <- temp.data
  }
  
  
  # Subset the clinical data so that it corresponds to individuals
  # in the intersected data
  intersectedIDs <- colnames(intersectedData[[1]])
  intersectedIDs <- transformIDFormat(intersectedIDs)
  # Find the position of the patient IDs in the clinical data
  clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
  clinicalInd <- match(intersectedIDs, clinicalIDs)
  clinicalData <- clinicalData[clinicalInd, ]
  
  # Subset clinical data and completeData by clinicalData
  clinicalData <- clinicalData[rowSums(is.na(clinicalData)) < ncol(clinicalData),]
  clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
  intersectedIndex <- match(clinicalIDs, intersectedIDs)
  
  clinicalData$days_to_death <- as.numeric(clinicalData$days_to_death)
  clinicalData$days_to_last_followup <- as.numeric(clinicalData$days_to_last_followup)
  ###### Run combat function 
  for (i in 1:numViews) {
    temp.intersectedData <- intersectedData[[i]]
    temp.intersectedData <- temp.intersectedData[,intersectedIndex]
    temp.2.intersectedData <- temp.intersectedData[rowSums(temp.intersectedData, na.rm = T) != 0,]
    temp.modcombat <- model.matrix(~1, data = clinicalData)
    temp.batch <- clinicalData$gender
    temp_combat = ComBat(dat=temp.2.intersectedData, batch=temp.batch, mod=temp.modcombat, par.prior=TRUE, prior.plots=FALSE)
    intersectedData[[i]] <- temp_combat
  }
  
  
} else {
  
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
  
  # Extract all cases which appear in at least one of the data types
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
  
  # Extract all cases which appear in all of the data types
  intersectedData <- columnIntersection(cases)
  
  # Subset the clinical data so that it corresponds to individuals
  # in the intersected data
  
  for (i in 1:3) {
    temp.data  <- intersectedData[[i]]
    temp.names <- transformUnionId(colnames(temp.data))
    colnames(temp.data) <- temp.names
    temp.data <- temp.data[, !duplicated(colnames(temp.data))]
    intersectedData[[i]] <- temp.data
  }
  
  # Subset the clinical data so that it corresponds to individuals
  # in the intersected data
  intersectedIDs <- colnames(intersectedData[[1]])
  intersectedIDs <- transformIDFormat(intersectedIDs)
  # Find the position of the patient IDs in the clinical data
  clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
  clinicalInd <- match(intersectedIDs, clinicalIDs)
  clinicalData <- clinicalData[clinicalInd, ]
  intersectedIndex <- match(intersectedIDs, intersectedIDs)
  
}
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

unionStat <- rowStatistics(unionData)
unionData <- normalizeData(unionData, unionStat)

intersectedStat <- rowStatistics(intersectedData)
intersectedData <- normalizeData(intersectedData, intersectedStat)
######################################################################
# Evaluate the performance of the similarity methods

sampleRows <- FALSE
similarityMethods <- c(selfSimilarity, medianSimilarity,
  regressionSimilarity)

similarityData <- list(unionData, intersectedData)

# Read in the cluster labels for the complete data
fileName <- paste(cancerInd, "_3.txt", sep="")
filePath <- paste(testFolder, "../cluster_complete_data",
                  "Results/Labels", fileName, sep="/")
completeLabelsInt <- scan(filePath)
#completeLabelsInt <- completeLabelsInt[intersectedIndex]


# Extend the completeLabels
nSamples <- ncol(unionData[[1]])
repTimes <- ceiling(nSamples/length(completeLabelsInt))
completeLabelsFull <- rep.int(completeLabelsInt, repTimes)[1:nSamples]

# Save the results of clustering the union and intersected data
for (i in 1:length(similarityMethods)) {
  similarity <- similarityMethods[[i]]
  
  for (j in 1:length(similarityData)) {
    data <- similarityData[[j]]
    tcgaID <- function(x) colnames(x[[1]])
    dataInd <- match(tcgaID(data), tcgaID(unionData))
    intersectInd <- tcgaID(data) %in% tcgaID(intersectedData)
    similarityResults <- evaluateSimilarityInt(similarity, data,
                                               completeLabelsFull[dataInd],
                                               clinicalData[dataInd, ],
                                               intersectInd,
                                               sampleRows)
    writeResults(c(argv, i, j, similarityResults), similarityFile)
  }
}
