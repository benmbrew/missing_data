############ This script is used to compare the summary stats in each cluster of the complete data
# to the summary stats of the patients left out of the column intersection. 
# Steps:
# 1) Find patients left out of columnIntersection
# 2) get summary stats of these patients 
# 3) get summary stats of complete data for each cluster. Both with imputation and not
# 4) compare clusters summary stats to left out summary stats. 

# Look at particular case: KIRC, Randomimputation and SNF clustering. 


######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)
library(dplyr)

######################################################################
# Initialize folders and load data 
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "evaluate_original_imputation", sep="/")
resultsFolder <- paste(testFolder, "Results", sep="/")
numFeat <- 2000
cancer <- 'KIRC'
cancerInd <- 'KIRC'

dataTypes <- c("methyl", "mirna", "mrna")
source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

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

###################################################################

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

# Union data
unionData <- columnUnion(cases)
unionIDs <- colnames(unionData[[1]])
unionIDs <- transformIDFormat(unionIDs)
clinicalIDsUnion <- as.character(clinicalData$bcr_patient_barcode)
clinicalIndUnion <- match(unionIDs, clinicalIDsUnion)
clinicalDataUnion <- clinicalData[clinicalIndUnion, ]

# Complete data 
completeData <- columnIntersection(cases)
completeIDs <- colnames(completeData[[1]])
completeIDs <- transformIDFormat(completeIDs)
clinicalIDsComplete <- as.character(clinicalData$bcr_patient_barcode)
clinicalIndComplete <- match(completeIDs, clinicalIDsComplete)
clinicalDataComplete <- clinicalData[clinicalIndComplete, ]

# get patients left our of the column intersection, that is patients in union, not in intersection.
clinicalExtra <- clinicalDataUnion[!clinicalDataUnion$bcr_patient_barcode %in% clinicalDataComplete$bcr_patient_barcode,]

# Summary stats of this group 
statsExtra <- summary(clinicalExtra)

################################################################################################
# Cluster the complete data with SNF 5 clusters, not random removal or imptutation 

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

completeStat <- rowStatistics(completeData)
completeData <- normalizeData(completeData, completeStat)

unionStat <- rowStatistics(unionData)
unionData <- normalizeData(unionData, unionStat)
###### Cluster with SNF
sampleRows = FALSE
numClus = 5
labelsComplete <- SNFClustering(completeData, numClus, sampleRows)

# merge lables with ids
labelsComplete <- as.data.frame(cbind(labelsComplete, transformIDFormat(completeIDs)))
colnames(labelsComplete)[2] <- 'id'
colnames(clinicalDataComplete)[1] <- 'id'
completeClinicalLabels <- left_join(labelsComplete, clinicalDataComplete, by = 'id')

# Subset by labels and take summary stats 
complete1 <- completeClinicalLabels[completeClinicalLabels$labelsComplete == 1,]
complete2 <- completeClinicalLabels[completeClinicalLabels$labelsComplete == 2,]
complete3 <- completeClinicalLabels[completeClinicalLabels$labelsComplete == 3,]
complete4 <- completeClinicalLabels[completeClinicalLabels$labelsComplete == 4,]
complete5 <- completeClinicalLabels[completeClinicalLabels$labelsComplete == 5,]

completeStats1 <- summary(complete1)
completeStats2 <- summary(complete2)
completeStats3 <- summary(complete3)
completeStats4 <- summary(complete4)
completeStats5 <- summary(complete5)

##################################################################
# randomly remove from complete data using generateincompleteintersection
# function. cluster these with SNF and compare clusters to statsExtra

incompleteIntersection <- generateIncompleteIntersection(cases)
incompleteData <- incompleteIntersection$incomplete

# impute using randomimputation and cluster that using SNF 
imputedData <- randomImputation(incompleteData, sampleRows)

# Use SNF on imputedData
sampleRows = FALSE
numClus = 5
labelsImpute <- SNFClustering(imputedData, numClus, sampleRows)

# merge lables with ids
labelsImpute <- as.data.frame(cbind(labelsImpute, transformIDFormat(completeIDs)))
colnames(labelsImpute)[2] <- 'id'
imputeClinicalLabels <- left_join(labelsImpute, clinicalDataComplete, by = 'id')

# Subset by labels and take summary stats 
completeImpute1 <- imputeClinicalLabels[imputeClinicalLabels$labelsImpute == 1,]
completeImpute2 <- imputeClinicalLabels[imputeClinicalLabels$labelsImpute == 2,]
completeImpute3 <- imputeClinicalLabels[imputeClinicalLabels$labelsImpute == 3,]
completeImpute4 <- imputeClinicalLabels[imputeClinicalLabels$labelsImpute == 4,]
completeImpute5 <- imputeClinicalLabels[imputeClinicalLabels$labelsImpute == 5,]

completeImputeStats1 <- summary(completeImpute1)
completeImputeStats2 <- summary(completeImpute2)
completeImputeStats3 <- summary(completeImpute3)
completeImputeStats4 <- summary(completeImpute4)
completeImputeStats5 <- summary(completeImpute5)

###############################################

# impute using randomimputation and cluster that using SNF 
imputedUnion <- randomImputation(unionData, sampleRows)

# Use SNF on imputedData
sampleRows = FALSE
numClus = 5
labelsImputedUnion <- SNFClustering(imputedUnion, numClus, sampleRows)

# merge lables with ids
labelsImputedUnion <- as.data.frame(cbind(labelsImputedUnion, transformIDFormat(unionIDs)))
colnames(labelsImputedUnion)[2] <- 'id'
colnames(clinicalData)[1] <- 'id'
imputedUnionClinicalLabels <- left_join(labelsImputedUnion, clinicalData, by = 'id')

# Subset by labels and take summary stats 
union1 <- imputedUnionClinicalLabels[imputedUnionClinicalLabels$labelsImputedUnion == 1,]
union2 <- imputedUnionClinicalLabels[imputedUnionClinicalLabels$labelsImputedUnion == 2,]
union3 <- imputedUnionClinicalLabels[imputedUnionClinicalLabels$labelsImputedUnion == 3,]
union4 <- imputedUnionClinicalLabels[imputedUnionClinicalLabels$labelsImputedUnion == 4,]
union5 <- imputedUnionClinicalLabels[imputedUnionClinicalLabels$labelsImputedUnion == 5,]

unionStats1 <- summary(union1)
unionStats2 <- summary(union2)
unionStats3 <- summary(union3)
unionStats4 <- summary(union4)
unionStats5 <- summary(union5)

##################################################################
# Look at id groups in complete data (random imputation) and compare to id groups 
# in original data (random) imputation) (completeImpute and union)

# try to merge clinical imputed and clinical union, keeping just the intersection
# first subset imputedUnionClinicalLabels by the intersection in IDs of imputeClinicalLabels

intersectLabels <- imputedUnionClinicalLabels[imputedUnionClinicalLabels$id %in% 
                                                imputeClinicalLabels$id,]

# now left_join intersectLabels and imputeClinicalLabels
intersectLabels <- left_join(intersectLabels, imputeClinicalLabels)

# Keep only labels and ids
intersectLabels <- intersectLabels[, c("id", "labelsImputedUnion", "labelsImpute")]

common_clusters <- matrix(,0,5)
cluster_compare <- matrix(,0,5)
for(i in 1:numClus){
  sub_union <- intersectLabels[intersectLabels$labelsImputedUnion == i,]
  for(j in 1:numClus){
    sub_intersect <- intersectLabels[intersectLabels$labelsImpute == j,]
    cluster_compare[j] <- sum(sub_union$id %in% sub_intersect$id)
  }
  common_clusters <- rbind(common_clusters, cluster_compare)
}


