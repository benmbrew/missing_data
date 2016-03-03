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
#   
#   ####################################################################################
# #   Normalize the features in the data sets.
# #   Normalization is performed before imputation and we expect that the
# #   data will still be normalized after imputation (before clustering).
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
groupByVital <- function(data, title) {
  
  temp <- data %>%
    group_by(vital_status) %>%
    summarise(counts = n())
  
  temp$percent <- (temp$counts/nrow(data))*100
  
  ggplot(data = temp, aes(vital_status, percent)) + geom_bar(stat = 'identity') + ggtitle(title)
  
}

groupByVital(clinicalDataCompleteBRCA, title = 'BRCA Complete')
groupByVital(clinBrca, title = 'BRCA Not in Intersection')
groupByVital(clinicalDataCompleteKIRC, title = 'KIRC Complete')
groupByVital(clinKirc, title = 'KIRC Not in Intersection')
groupByVital(clinicalDataCompleteLIHC, title = 'LIHC Complete')
groupByVital(clinLihc, title = 'LIHC Not in Intersection')
groupByVital(clinicalDataCompleteLUAD, title = 'LUAD Complete')
groupByVital(clinLuad, title = 'LUAD Not in Intersection')


# group by gender and return a bargraph of alive, dead, and NAs percentage for complete and left out clinical ids.
groupByGender <- function(data, title) {
  
  temp <- data %>%
    group_by(gender) %>%
    summarise(counts = n())
  
  temp$percent <- (temp$counts/nrow(data))*100
  
  ggplot(data = temp, aes(gender, percent)) + geom_bar(stat = 'identity') + ggtitle(title)
  
}

groupByGender(clinicalDataCompleteBRCA, title = 'BRCA Complete')
groupByGender(clinBrca, title = 'BRCA Not in Intersection')
groupByGender(clinicalDataCompleteKIRC, title = 'KIRC Complete')
groupByGender(clinKirc, title = 'KIRC Not in Intersection')
groupByGender(clinicalDataCompleteLIHC, title = 'LIHC Complete')
groupByGender(clinLihc, title = 'LIHC Not in IntersectionUnion')
groupByGender(clinicalDataCompleteLUAD, title = 'LUAD Union')
groupByGender(clinLuad, title = 'LUAD Not in Intersection')

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
daysToDeath(clinBrca, title = 'BRCA Not in Intersection', column = 'Days To Death')

# KIRC
daysToDeath(clinicalDataCompleteKIRC, title = 'KIRC Complete', column = 'Days To Death')
daysToDeath(clinKirc, title = 'KIRC Not in Intersection', column = 'Days To Death')

# LIHC
daysToDeath(clinicalDataCompleteLIHC, title = 'LIHC Complete', column = 'Days To Death')
daysToDeath(clinLihc, title = 'LIHC Not in Intersection', column = 'Days To Death')

#LUAD
daysToDeath(clinicalDataCompleteLUAD, title = 'LUAD Complete', column = 'Days To Death')
daysToDeath(clinLuad, title = 'LUAD Not in Intersection', column = 'Days To Death')


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
daysToLastFollowUp(clinBrca, title = 'BRCA Not in Intersection', column = 'Days To Last Follow Up')

# KIRC
daysToLastFollowUp(clinicalDataCompleteKIRC, title = 'KIRC Complete', column = 'Days To Last Follow Up')
daysToLastFollowUp(clinKirc, title = 'KIRC Not in Intersection', column = 'Days To Last Follow Up')

# LIHC
daysToLastFollowUp(clinicalDataCompleteLIHC, title = 'LIHC Complete', column = 'Days To Last Follow Up')
daysToLastFollowUp(clinLihc, title = 'LIHC Not in Intersection', column = 'Days To Last Follow Up')

#LUAD
daysToLastFollowUp(clinicalDataCompleteLUAD, title = 'LUAD Complete', column = 'Days To Last Follow Up')
daysToLastFollowUp(clinLuad, title = 'LUAD Not in Intersection', column = 'Days To Last Follow Up')


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

survObject(clinicalDataCompleteBRCA, 'BRCA Complete')
survObject(clinBrca, 'BRCA Not in Intersection')

survObject(clinicalDataCompleteKIRC, 'KIRC Complete')
survObject(clinKirc, 'KIRC Not in Intersection')

survObject(clinicalDataCompleteLIHC, 'LIHC')
survObject(clinLihc, 'LIHC Not in Intersection')

survObject(clinicalDataCompleteLUAD, 'LUAD')
survObject(clinLuad, 'LUAD Not in Intersection')
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
  if(column == 'vital_status'){
    
    temp.diff_BRCA <- round((nrow(clinicalDataCompleteBRCA[which(clinicalDataCompleteBRCA$vital_status == 'alive'),])/
                              nrow(clinicalDataCompleteBRCA))*100 - 
                              (nrow(clinBrca[which(clinBrca$vital_status == 'alive'),])/nrow(clinBrca))*100, 2)
    temp.diff_KIRC <- round((nrow(clinicalDataCompleteKIRC[which(clinicalDataCompleteKIRC$vital_status == 'alive'),])/
                              nrow(clinicalDataCompleteKIRC))*100 - 
                              (nrow(clinKirc[which(clinKirc$vital_status == 'alive'),])/nrow(clinKirc))*100, 2)
    temp.diff_LIHC <- round((nrow(clinicalDataCompleteLIHC[which(clinicalDataCompleteLIHC$vital_status == 'alive'),])/
                               nrow(clinicalDataCompleteLIHC))*100 - 
                              (nrow(clinLihc[which(clinLihc$vital_status == 'alive'),])/nrow(clinLihc))*100, 2)
    temp.diff_LUAD <- round((nrow(clinicalDataCompleteLUAD[which(clinicalDataCompleteLUAD$vital_status == 'alive'),])/
                               nrow(clinicalDataCompleteLUAD))*100 - 
                              (nrow(clinLuad[which(clinLuad$vital_status == 'alive'),])/nrow(clinLuad))*100, 2) 
  
  } else if (column == 'gender') {
    temp.diff_BRCA <- round((nrow(clinicalDataCompleteBRCA[which(clinicalDataCompleteBRCA$gender == 'female'),])/
                               nrow(clinicalDataCompleteBRCA))*100 - 
                              (nrow(clinBrca[which(clinBrca$gender == 'female'),])/nrow(clinBrca))*100, 2)
    temp.diff_KIRC <- round((nrow(clinicalDataCompleteKIRC[which(clinicalDataCompleteKIRC$gender == 'female'),])/
                               nrow(clinicalDataCompleteKIRC))*100 - 
                              (nrow(clinKirc[which(clinKirc$gender == 'female'),])/nrow(clinKirc))*100, 2)
    temp.diff_LIHC <- round((nrow(clinicalDataCompleteLIHC[which(clinicalDataCompleteLIHC$gender == 'female'),])/
                               nrow(clinicalDataCompleteLIHC))*100 - 
                              (nrow(clinLihc[which(clinLihc$gender == 'female'),])/nrow(clinLihc))*100, 2)
    temp.diff_LUAD <- round((nrow(clinicalDataCompleteLUAD[which(clinicalDataCompleteLUAD$gender == 'female'),])/
                               nrow(clinicalDataCompleteLUAD))*100 - 
                              (nrow(clinLuad[which(clinLuad$gender == 'female'),])/nrow(clinLuad))*100, 2)
  
  } else {
  temp.diff_BRCA <- mean(clinicalDataCompleteBRCA[, column], na.rm= T) - mean(clinBrca[, column], na.rm= T)
  temp.diff_KIRC <- mean(clinicalDataCompleteKIRC[, column], na.rm= T) - mean(clinKirc[, column], na.rm= T)
  temp.diff_LIHC <- mean(clinicalDataCompleteLIHC[, column], na.rm= T) - mean(clinLihc[, column], na.rm= T)
  temp.diff_LUAD <- mean(clinicalDataCompleteLUAD[, column], na.rm= T) - mean(clinLuad[, column], na.rm= T)
  }
  
  temp <- data.frame(x = c('BRCA', 'KIRC', 'LIHC', 'LUAD'), y = c(temp.diff_BRCA, temp.diff_KIRC, temp.diff_LIHC, temp.diff_LUAD))
  
  ggplot(temp, aes(x, y)) + geom_bar(stat = 'identity') + ggtitle(title)

}

Difference(column = 'days_to_death', title = 'Difference in days_to_death between complete and not complete')
Difference(column = 'days_to_last_followup', title = 'Difference in days_to_last_followup between complete and not complete')
Difference(column = 'survTime', title = 'Difference in Survial Time between complete and not complete')
Difference(column = 'vital_status', title = 'Difference in Percent Alive between complete and not complete')
Difference(column = 'gender', title = 'Difference in Percent Female between complete and not complete')


#############################################################################################
# # Run PCA on complete data without combat and with combat on KIRC
# library(sva)
# 
# # transform patient IDs to the clinical ID format 
# transformIDFormat <- function(x){
#   x <- substr(x, 1, 12)
#   x <- gsub('.', '-', x, fixed = TRUE)
#   x <- tolower(x)
#   
#   return(x)
# }
# 
# runCombat <- function(completeData, clinicalData, numViews = 3) {
#   
#   for (i in 1:numViews) {
#     # Subset clinical data and completeData by clinicalData
#     completeIDs <- colnames(completeData[[1]])
#     completeIDs <- transformIDFormat(completeIDs)
#     clinicalData <- clinicalData[rowSums(is.na(clinicalData)) < ncol(clinicalData),]
#     clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
#     completeInd <- match(clinicalIDs, completeIDs)
#     temp.completeData <- completeData[[i]]
#     temp.completeData <- temp.completeData[,completeInd]
#     temp.2.completeData <- temp.completeData[rowSums(temp.completeData) != 0,]
#     temp.modcombat <- model.matrix(~1, data = clinicalData)
#     temp.batch <- clinicalData$gender
#     temp_combat = ComBat(dat=temp.2.completeData, batch=temp.batch, mod=temp.modcombat, par.prior=TRUE, prior.plots=FALSE)
#     completeData[[i]] <- temp_combat
#   }
#   return(completeData)
# }
# 
# dataMerge <- function(data, clinical){
#   data <- as.data.frame(t(data))
#   data <- cbind('id' = row.names(data), data)
#   data$id <- transformIDFormat(data$id)
#   row.names(data) <- NULL
#   names(clinical)[2] <- 'id'
#   data <- left_join(data, clinical, by = 'id')
#   data$days_to_death[is.na(as.numeric(data$days_to_death))] <- max(clinical$days_to_death, na.rm = T)
#   return(data)
# }
# 
# pca <- function(data){
#   data <- data[!is.na(data$days_to_death),]
#   data <- data[data$days_to_death > 0,]
#   data_length <- (ncol(data)-8) 
#   pca <- prcomp(data[,2:data_length])
#   return(pca)
# }
# 
# # get gene info
# KIRCCompleteGene <- KIRCComplete[[1]]
# 
# # run combat on KircCompleteGene
# KircCombat <- runCombat(KIRCCompleteGene, clinicalDataCompleteKIRC)
# 
# # get methylation from gene info and combat
# KircMethyl <- KIRCCompleteGene[[1]]
# KircMethylCombat <- KircCombat[[1]]
# 
# # Merge both with clinical data
# KircMethylFull <- dataMerge(KircMethyl, clinicalDataCompleteKIRC)
# KircMethylCombatFull <- dataMerge(KircMethylCombat, clinicalDataCompleteKIRC)
# 
# # run pca on both data sets 
# KircMethylFullPca <- pca(KircMethylFull)
# KircMethylCombatPca <- pca(KircMethylCombatFull)
# 
# # Plot PCA
# pcaPlot <- function(pca, name){
# #   min <- min(min(pca$x[,1]), pca$x[,2])
# #   max <- max(max(pca$x[,1]), pca$x[,2])
#   plot(pca$x[,1], 
#        pca$x[,2],
#        xlab = 'PCA 1',
#        ylab = 'PCA 2',
#        cex = 1,
#        main = name,
#        pch = 16
# #        xlim= c(min, max),
# #        ylim = c(min, max)
#   )
#   abline(v = c(0,0),
#          h = c(0,0))
#   
# }
# 
# pcaPlot(KircMethylFullPca, name = 'PCA')
# pcaPlot(KircMethylCombatPca, name = 'PCA')
# 
# #########################################################################################
