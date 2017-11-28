#argv <- as.numeric(commandArgs(T))

######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)

######################################################################
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")

testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features", 
                    "cluster_complete_data", sep="/")
resultsFolder <- paste(testFolder, "Results", sep="/")

dataTypes <- c("methyl", "mirna", "mrna")
source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

### Load in complete labels of all cancers. 
# Need to compare this to the union labels

cancer <- c("BRCA", "KIRC", "LIHC", "LUAD")
dataTypes <- c("methyl", "mirna", "mrna")

clusteringMethods <- c(hierarchicalClustering, iClusterClustering,
                       SNFClustering)
BRCA <- vector('list', 3)
KIRC <- vector('list', 3)
LIHC <- vector('list', 3)
LUAD <- vector('list', 3)

setwd('/home/benbrew/Documents/NM_2015_ben/Scripts/06_two_thousand_features/cluster_complete_data/Results/Labels')
for(i in 1:3){
  BRCA[[i]] <- as.matrix(read.table(paste('1_' ,i, '.txt', sep = '')))
}
for(i in 1:3){
  KIRC[[i]] <- as.matrix(read.table(paste('2_' ,i, '.txt', sep = '')))
}
for(i in 1:3){
  LIHC[[i]] <- as.matrix(read.table(paste('3_' ,i, '.txt', sep = '')))
}
for(i in 1:3){
  LUAD[[i]] <- as.matrix(read.table(paste('4_' ,i, '.txt', sep = '')))
}

BRCA <- do.call('rbind', BRCA)
KIRC <- do.call('rbind', KIRC)
LIHC <- do.call('rbind', LIHC)
LUAD <- do.call('rbind', LUAD)

BRCA <- as.data.frame(t(BRCA))
KIRC <- as.data.frame(t(KIRC))
LIHC <- as.data.frame(t(LIHC))
LUAD <- as.data.frame(t(LUAD))

BRCA <- apply(BRCA, 2, as.factor)
KIRC <- apply(KIRC, 2, as.factor)
LIHC <- apply(LIHC, 2, as.factor)
LUAD <- apply(LUAD, 2, as.factor)

### examine KIRC ############################################################################
# first look at clinical 

# compare intersected clin and union clin
load_cancer <- function(cancer){
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

clinicalData <- loadData("clin")

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

completeData <- columnIntersection(cases)

# Subset the clinical data so that it corresponds to individuals
# in the complete data
completeIDs <- colnames(completeData[[1]])
completeIDs <- transformIDFormat(completeIDs)
# Find the position of the patient IDs in the clinical data
clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
clinicalInd <- match(completeIDs, clinicalIDs)
clinicalData <- clinicalData[clinicalInd, ]
return(clinicalData)
}
BRCA_clin_inter <- load_cancer(cancer = 'BRCA')
KIRC_clin_inter <- load_cancer(cancer = 'KIRC')
LIHC_clin_inter <- load_cancer(cancer = 'LIHC')
LUAD_clin_inter <- load_cancer(cancer = 'LUAD')

# the union 
setwd('/home/benbrew/Documents/NM_2015_ben/Data')
BRCA_clin <- read.table('BRCA_clin.txt', header = TRUE)
KIRC_clin <- read.table('KIRC_clin.txt', header = TRUE)
LIHC_clin <- read.table('LIHC_clin.txt', header = TRUE)
LUAD_clin <- read.table('LUAD_clin.txt', header = TRUE)

# biggest percent reduction
(nrow(BRCA_clin_inter)/nrow(BRCA_clin))*100
(nrow(KIRC_clin_inter)/nrow(KIRC_clin))*100
(nrow(LIHC_clin_inter)/nrow(LIHC_clin))*100
(nrow(LUAD_clin_inter)/nrow(LUAD_clin))*100
# KIRC has the biggest reduction (40%) from union to intersection
## Look at vital status, days to death, days to last followup, days to last known alive

#################vital status
summary(BRCA_clin$vital_status)
summary(BRCA_clin_inter$vital_status)
summary(KIRC_clin$vital_status)
summary(KIRC_clin_inter$vital_status)
summary(LIHC_clin$vital_status)
summary(LIHC_clin_inter$vital_status)
summary(LUAD_clin$vital_status)
summary(LUAD_clin_inter$vital_status)

# BRCA
nrow(BRCA_clin[BRCA_clin$vital_status == 'dead',])/
  nrow(BRCA_clin[BRCA_clin$vital_status == 'alive',])
nrow(BRCA_clin_inter[BRCA_clin_inter$vital_status == 'dead',])/
  nrow(BRCA_clin_inter[BRCA_clin_inter$vital_status == 'alive',])

# KIRC
nrow(KIRC_clin[KIRC_clin$vital_status == 'dead',])/
  nrow(KIRC_clin[KIRC_clin$vital_status == 'alive',])
nrow(KIRC_clin_inter[KIRC_clin_inter$vital_status == 'dead',])/
  nrow(KIRC_clin_inter[KIRC_clin_inter$vital_status == 'alive',])

# LIHC
nrow(LIHC_clin[LIHC_clin$vital_status == 'dead',])/
  nrow(LIHC_clin[LIHC_clin$vital_status == 'alive',])
nrow(LIHC_clin_inter[LIHC_clin_inter$vital_status == 'dead',])/
  nrow(LIHC_clin_inter[LIHC_clin_inter$vital_status == 'alive',])

# LUAD
nrow(LUAD_clin[LUAD_clin$vital_status == 'dead',])/
  nrow(LUAD_clin[LUAD_clin$vital_status == 'alive',])
nrow(LUAD_clin_inter[LUAD_clin_inter$vital_status == 'dead',])/
  nrow(LUAD_clin_inter[LUAD_clin_inter$vital_status == 'alive',])

## In union data, KIRC has highest percentage of dead people by far.

############### Days to death 
summary(BRCA_clin$days_to_death)
summary(BRCA_clin_inter$days_to_death)
summary(KIRC_clin$days_to_death)
summary(KIRC_clin_inter$days_to_death)
summary(LIHC_clin$days_to_death)
summary(LIHC_clin_inter$days_to_death)
summary(LUAD_clin$days_to_death)
summary(LUAD_clin_inter$days_to_death)

############### Days to last followup 
summary(BRCA_clin$days_to_last_followup)
summary(BRCA_clin_inter$days_to_last_followup)
summary(KIRC_clin$days_to_last_followup)
summary(KIRC_clin_inter$days_to_last_followup)
summary(LIHC_clin$days_to_last_followup)
summary(LIHC_clin_inter$days_to_last_followup)
summary(LUAD_clin$days_to_last_followup)
summary(LUAD_clin_inter$days_to_last_followup)
## KIRC has highest days to last follow up in both union and intersection

############# replace NAs in day to death with days to last follow up
# Replace missing survival times with days to last follow up
clinical_days <- function(clinicalData){
survTime <- as.numeric(clinicalData$days_to_death)
missingSurvInd <- is.na(survTime)
lastFollowup <- clinicalData$days_to_last_followup
survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
return(survTime)
}
BRCA_time <- clinical_days(BRCA_clin)
BRCA_int_time <- clinical_days(BRCA_clin_inter)
KIRC_time <- clinical_days(KIRC_clin)
KIRC_int_time <- clinical_days(KIRC_clin_inter)
LIHC_time <- clinical_days(LIHC_clin)
LIHC_int_time <- clinical_days(LIHC_clin_inter)
LUAD_time <- clinical_days(LUAD_clin)
LUAD_int_time <- clinical_days(LUAD_clin_inter)

summary(BRCA_time)
summary(BRCA_int_time)
summary(KIRC_time)
summary(KIRC_int_time)
summary(LIHC_time)
summary(LIHC_int_time)
summary(LUAD_time)
summary(LUAD_int_time)
# KIRC has highest time in both union and intersection

## Look at summary of data. 



######################################################################
# Initialize fixed variables

cancerTypes <- c("KIRC")
dataTypes <- c("methyl", "mirna", "mrna")
cancer <- 'KIRC'
source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

######################################################################
# Load the original data

cancer_cases <- function(cancer){

  numFeat <- 2000
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

completeData <- columnIntersection(cases)



# Extract all cases which appear in all of the data types in the
# incomplete data


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
return(completeData)
}

BRCA_complete <- cancer_cases(cancer = 'BRCA')
KIRC_complete <- cancer_cases(cancer = 'KIRC')
LIHC_complete <- cancer_cases(cancer = 'LIHC')
LUAD_complete <- cancer_cases(cancer = 'LUAD')

# BRCA_methyl <- as.data.frame(t(BRCA_complete[[1]]))
# BRCA_mirna <- as.data.frame(t(BRCA_complete[[2]]))
# BRCA_mrna <- as.data.frame(t(BRCA_complete[[3]]))
# 
# KIRC_methyl <- as.data.frame(t(KIRC_complete[[1]]))
# KIRC_mirna <- as.data.frame(t(KIRC_complete[[2]]))
# KIRC_mrna <- as.data.frame(t(KIRC_complete[[3]]))
# 
# LIHC_methyl <- as.data.frame(t(LIHC_complete[[1]]))
# LIHC_mirna <- as.data.frame(t(LIHC_complete[[2]]))
# LIHC_mrna <- as.data.frame(t(LIHC_complete[[3]]))
# 
# LUAD_methyl <- as.data.frame(t(LUAD_complete[[1]]))
# LUAD_mirna <- as.data.frame(t(LUAD_complete[[2]]))
# LUAD_mrna <- as.data.frame(t(LUAD_complete[[3]]))

####################################################
# Look at union data
# compare intersected clin and union clin
load_cases <- function(cancer){
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
  # Extract all cases which appear in at least one of the data types
  return(cases)
}  

BRCA_cases <- load_cases(cancer = 'BRCA')
KIRC_cases <- load_cases(cancer = 'KIRC')
LIHC_cases <- load_cases(cancer = 'LIHC')
LUAD_cases <- load_cases(cancer = 'LUAD')

membershipMatrix <- function(membershipList) {
  allMembers <- unique(unlist(membershipList), use.names=FALSE) # gets all ids across list of 3 data sets.
  sapply(membershipList, function(x) allMembers %in% x) 
}

BRCA_member <- membershipMatrix(lapply(BRCA_cases, colnames)) # Creates index of all patients and 
KIRC_member <- membershipMatrix(lapply(KIRC_cases, colnames)) # Creates index of all patients and 
LIHC_member <- membershipMatrix(lapply(LIHC_cases, colnames)) # Creates index of all patients and 
LUAD_member <- membershipMatrix(lapply(LUAD_cases, colnames)) # Creates index of all patients and 


