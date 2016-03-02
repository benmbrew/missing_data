###############################################################################
# This script will look at the difference between clinical data in the intersection and union.

# Load libraries
library(dplyr)
library(ggplot2)
###############################################################################
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Combat",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(testFolder, "Results", sep="/")

###############################################################################
# Initialize fixed variables
jvmGBLimit <- 8
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 12
numFeat <- 2000

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

################################################################################
# Load clinical data
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


loadClinData <- function(cancer) {
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
clinicalDataBRCA <- loadClinData(cancer = 'BRCA')
clinicalDataBRCA <- extractRelevantColumns(clinicalDataBRCA)
clinicalDataKIRC <- loadClinData(cancer = 'KIRC')
clinicalDataKIRC <- extractRelevantColumns(clinicalDataKIRC)
clinicalDataLIHC <- loadClinData(cancer = 'LIHC')
clinicalDataLIHC <- extractRelevantColumns(clinicalDataLIHC)
clinicalDataLUAD <- loadClinData(cancer = 'LUAD')
clinicalDataLUAD <- extractRelevantColumns(clinicalDataLUAD)

# remove patient from column names
removePatients <- function(data) {
  
  features <- colnames(data)
  split <- strsplit(features, '.', fixed = TRUE)
  keepSplit <- lapply(split, function(x) x[length(x)])
  features <- unlist(keepSplit)
  colnames(data) <- features
  return(data)
  
}

# apply the function
clinicalDataBRCA <- removePatients(clinicalDataBRCA)
clinicalDataKIRC <- removePatients(clinicalDataKIRC)
clinicalDataLIHC <- removePatients(clinicalDataLIHC)
clinicalDataLUAD <- removePatients(clinicalDataLUAD)


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
  
  ####################################################################################
#   Normalize the features in the data sets.
#   Normalization is performed before imputation and we expect that the
#   data will still be normalized after imputation (before clustering).
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


####################################################################################
# extract samples in union that are not in intersection for each data type 
clinBrca <- clinicalDataUnionBRCA[!clinicalDataUnionBRCA$bcr_patient_barcode %in% 
                                    clinicalDataCompleteBRCA$bcr_patient_barcode,]

clinKirc <- clinicalDataUnionKIRC[!clinicalDataUnionKIRC$bcr_patient_barcode %in% 
                                    clinicalDataCompleteKIRC$bcr_patient_barcode,]

clinLihc <- clinicalDataUnionLIHC[!clinicalDataUnionLIHC$bcr_patient_barcode %in% 
                                    clinicalDataCompleteLIHC$bcr_patient_barcode,]

clinLuad <- clinicalDataUnionLUAD[!clinicalDataUnionLUAD$bcr_patient_barcode %in% 
                                    clinicalDataCompleteLUAD$bcr_patient_barcode,]

####################################################################################
# summary statistics, histograms, for clinical data in the intersection

# change structure of data 
dataStructure <- function(data) {
  data$days_to_death <- as.numeric(as.character(data$days_to_death))
  data$days_to_last_followup <- as.numeric(as.character(data$days_to_last_followup))
  data$days_to_last_known_alive <- as.numeric(as.character(data$days_to_last_known_alive))
  return(data)
}

# Apply the function to the clinical data in the intersection and the clinical data not in the intersection. 
clinBrca <- dataStructure(clinBrca)
clinicalDataCompleteBRCA <- dataStructure(clinicalDataCompleteBRCA) 
clinKirc <- dataStructure(clinKirc)
clinicalDataCompleteKIRC <- dataStructure(clinicalDataCompleteKIRC) 
clinLihc <- dataStructure(clinLihc)
clinicalDataCompleteLIHC <- dataStructure(clinicalDataCompleteLIHC) 
clinLuad <- dataStructure(clinLuad)
clinicalDataCompleteLUAD <- dataStructure(clinicalDataCompleteLUAD) 


# dimension of data 
dim(clinicalDataUnionBRCA) #1102
dim(clinicalDataCompleteBRCA) # 716
dim(clinicalDataUnionKIRC) # 542
dim(clinicalDataCompleteKIRC) #315
dim(clinicalDataUnionLIHC) # 377
dim(clinicalDataCompleteLIHC) #284
dim(clinicalDataUnionLUAD) # 529
dim(clinicalDataCompleteLUAD) # 430

#########################################################################################
# Barplots of all vital_status and gender (with NAs) for complete data (clinicalDataComplete) and not complete (clin)


# group by vital status and return a bargraph of alive, dead, and NAs percentage for complete and left out clinical ids.
groupByVital <- function(data) {
  
  temp <- data %>%
    group_by(vital_status) %>%
    summarise(counts = n())
  
  temp$percent <- (temp$counts/nrow(data))*100
  
  ggplot(data = temp, aes(vital_status, percent)) + geom_bar(stat = 'identity')
  
}

groupByVital(clinicalDataCompleteBRCA)
groupByVital(clinBrca)
groupByVital(clinicalDataCompleteKIRC)
groupByVital(clinKirc)
groupByVital(clinicalDataCompleteLIHC)
groupByVital(clinLihc)
groupByVital(clinicalDataCompleteLUAD)
groupByVital(clinLuad)


# group by gender and return a bargraph of alive, dead, and NAs percentage for complete and left out clinical ids.
groupByGender <- function(data) {
  
  temp <- data %>%
    group_by(gender) %>%
    summarise(counts = n())
  
  temp$percent <- (temp$counts/nrow(data))*100
  
  ggplot(data = temp, aes(gender, percent)) + geom_bar(stat = 'identity')
  
}

groupByGender(clinicalDataCompleteBRCA)
groupByGender(clinBrca)
groupByGender(clinicalDataCompleteKIRC)
groupByGender(clinKirc)
groupByGender(clinicalDataCompleteLIHC)
groupByGender(clinLihc)
groupByGender(clinicalDataCompleteLUAD)
groupByGender(clinLuad)

##############################################################################################
# Histograms of days_to_death 
daysToDeath <- function(data, title, column) {
  hist(data$days_to_death, main = title, xlab = column)
  legend('topright', legend = c(paste0(round(length(which(is.na(data$days_to_death)))/
                                               nrow(data)*100), '% NA'),
                                paste0('Vital Status', ':'),
                                paste0(round((nrow(data[which(data$vital_status == 'alive'),])/nrow(data))*100), '% alive'),
                                paste0(round((nrow(data[which(data$vital_status == 'dead'),])/nrow(data))*100), '% dead')), bty = 'n') 
}

# BRCA
daysToDeath(clinicalDataCompleteBRCA, title = 'BRCA Complete', column = 'Days To Death')
daysToDeath(clinBrca, title = 'BRCA Union', column = 'Days To Death')

# KIRC
daysToDeath(clinicalDataCompleteKIRC, title = 'KIRC Complete', column = 'Days To Death')
daysToDeath(clinKirc, title = 'KIRC Union', column = 'Days To Death')

# LIHC
daysToDeath(clinicalDataCompleteLIHC, title = 'LIHC Complete', column = 'Days To Death')
daysToDeath(clinLihc, title = 'LIHC Union', column = 'Days To Death')

#LUAD
daysToDeath(clinicalDataCompleteLUAD, title = 'LUAD Complete', column = 'Days To Death')
daysToDeath(clinLuad, title = 'LUAD Union', column = 'Days To Death')


######################################################################################################
# Histograms of days_to_last_followup

daysToLastFollowUp <- function(data, title, column) {
  hist(data$days_to_last_followup, main = title, xlab = column)
  legend('topright', legend = c(paste0(round(length(which(is.na(data$days_to_last_followup)))/
                                               nrow(data)*100), '% NA'),
                                paste0('Vital Status', ':'),
                                paste0(round((nrow(data[which(data$vital_status == 'alive'),])/nrow(data))*100), '% alive'),
                                paste0(round((nrow(data[which(data$vital_status == 'dead'),])/nrow(data))*100), '% dead')), bty = 'n') 
}
# BRCA
daysToLastFollowUp(clinicalDataCompleteBRCA, title = 'BRCA Complete', column = 'Days To Last Follow Up')
daysToLastFollowUp(clinBrca, title = 'BRCA Union', column = 'Days To Last Follow Up')

# KIRC
daysToLastFollowUp(clinicalDataCompleteKIRC, title = 'KIRC Complete', column = 'Days To Last Follow Up')
daysToLastFollowUp(clinKirc, title = 'KIRC Union', column = 'Days To Last Follow Up')

# LIHC
daysToLastFollowUp(clinicalDataCompleteLIHC, title = 'LIHC Complete', column = 'Days To Last Follow Up')
daysToLastFollowUp(clinLihc, title = 'LIHC Union', column = 'Days To Last Follow Up')

#LUAD
daysToLastFollowUp(clinicalDataCompleteLUAD, title = 'LUAD Complete', column = 'Days To Last Follow Up')
daysToLastFollowUp(clinLuad, title = 'LUAD Union', column = 'Days To Last Follow Up')


##########################################################################################################
# Histogram of surv object
survObject <- function(data, title) {
  
  survTime <- data$days_to_death
  deathStatus <- data$vital_status == "dead"
  
  # Replace missing survival times with days to last follow up
  missingSurvInd <- is.na(survTime)
  lastFollowup <- data$days_to_last_followup
  survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
  hist(survTime, main = title, xlab = 'Survival')
  legend('topright', legend = c(paste0(round(length(which(is.na(survTime)))/
                                             length(survTime)*100), '% NA'), 
                                paste0(round((nrow(data[which(data$vital_status == 'alive'),])/nrow(data))*100), '% alive'),
                                paste0(round((nrow(data[which(data$vital_status == 'dead'),])/nrow(data))*100), '% dead')), bty = 'n')
  
  
}

survObject(clinicalDataCompleteBRCA, 'BRCA')
survObject(clinBrca, 'BRCA')

survObject(clinicalDataCompleteKIRC, 'KIRC')
survObject(clinKirc, 'KIRC')

survObject(clinicalDataCompleteLIHC, 'LIHC')
survObject(clinLihc, 'LIHC')

survObject(clinicalDataCompleteLUAD, 'LUAD')
survObject(clinLuad, 'LUAD')
###########################################################################################
# Create a function to add survTime to clinicalDataComplete.... and clin....
addSurv <- function(data) {
  data$survTime <- data$days_to_death
  missingSurvInd <- is.na(data$survTime)
  lastFollowup <- data$days_to_last_followup
  data$survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
  return(data)
}

# apply the function. now each data set has a survTime column
clinicalDataCompleteBRCA <- addSurv(clinicalDataCompleteBRCA)
clinBrca <- addSurv(clinBrca)
clinicalDataCompleteKIRC <- addSurv(clinicalDataCompleteKIRC)
clinKirc <- addSurv(clinKirc)
clinicalDataCompleteLIHC <- addSurv(clinicalDataCompleteLIHC)
clinLihc <- addSurv(clinLihc)
clinicalDataCompleteLUAD <- addSurv(clinicalDataCompleteLUAD)
clinLuad <- addSurv(clinLuad)



###########################################################################################
# Look at relative difference between percentiles. 
# percent difference for each 
Difference <- function(column, title) {
  
  temp.diff_BRCA <- mean(clinicalDataCompleteBRCA[, column], na.rm= T) - mean(clinBrca[, column], na.rm= T)
  temp.diff_KIRC <- mean(clinicalDataCompleteKIRC[, column], na.rm= T) - mean(clinKirc[, column], na.rm= T)
  temp.diff_LIHC <- mean(clinicalDataCompleteLIHC[, column], na.rm= T) - mean(clinLihc[, column], na.rm= T)
  temp.diff_LUAD <- mean(clinicalDataCompleteLUAD[, column], na.rm= T) - mean(clinLuad[, column], na.rm= T)
  
  temp <- data.frame(x = c('BRCA', 'KIRC', 'LIHC', 'LUAD'), y = c(temp.diff_BRCA, temp.diff_KIRC, temp.diff_LIHC, temp.diff_LUAD))
  
  ggplot(temp, aes(x, y)) + geom_bar(stat = 'identity') + ggtitle(title)

}

Difference(column = 'days_to_death', title = 'Difference in days_to_death between complete and union')
Difference(column = 'days_to_last_followup', title = 'Difference in days_to_last_followup between complete and union')
Difference(column = 'survTime', title = 'Difference in Survial Time between complete and union')
