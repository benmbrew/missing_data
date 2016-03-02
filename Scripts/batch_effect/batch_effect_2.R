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
  
  #   featureSubsetIndices <- function(cases, subsetSize=numFeat) {
  #     numViews <- length(cases)
  #     featureSubsetInd <- vector("list", numViews)
  #     
  #     for (v in 1:numViews) {
  #       # Calculate the t-test p-value for each feature, grouped by cases
  #       # and controls
  #       numFeatures <- nrow(cases[[v]])
  #       pval <- sapply(1:numFeatures,
  #                      function(i) t.test(cases[[v]][i, ],
  #                                         controls[[v]][i, ])$p.value)
  #       
  #       # Subset the data keeping the features with the smallest p-values
  #       ind <- order(pval)
  #       featureSubsetInd[[v]] <- ind[1:min(subsetSize, numFeatures)]
  #     }
  #     
  #     return(featureSubsetInd)
  #   }
  #   
  #   subsetData <- function(data, ind) {
  #     
  #     for (v in 1:length(data)) {
  #       data[[v]] <- data[[v]][ind[[v]], ]
  #     }
  #     
  #     return(data)
  #   }
  #   
  #  if (complete){
  #     completeInd <- featureSubsetIndices(completeData)
  #     completeData <- subsetData(completeData, completeInd)
  #     } else {
  #     unionInd <- featureSubsetIndices(unionData)
  #     unionData <- subsetData(unionData, unionInd)
  #     }
  # #   
  #   ####################################################################################
  # #   Normalize the features in the data sets.
  # #   Normalization is performed before imputation and we expect that the
  # #   data will still be normalized after imputation (before clustering).
  #     rowStatistics <- function(cases){
  #       num_views <- length(cases)
  #       row_stats <- vector('list', num_views)
  #       
  #       for(v in 1:num_views){
  #         #calculate the row means and std deviations 
  #         row_mean <- apply(cases[[v]], 1, mean, na.rm = T)
  #         row_sd <- apply(cases[[v]], 1, sd, na.rm = T)
  #         constant_ind <- row_sd == 0
  #         row_sd[constant_ind] <- 1
  #         row_stats[[v]] <- list(mean = row_mean, sd = row_sd, ind = constant_ind)
  #       }
  #       return(row_stats)
  #     }
  #     
  #     normalizeData <- function(data, stat){
  #       for(v in 1:length(data)) {
  #         data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
  #         data[[v]] <- data[[v]][!stat[[v]]$ind, ]
  #       }
  #       return(data)
  #     }
  #     
  #    if (complete) {
  #       completeStat <- rowStatistics(completeData)
  #       completeData <- normalizeData(completeData, completeStat)
  #     } else {
  #     unionStat <- rowStatistics(unionData)
  #     unionData <- normalizeData(unionData, unionStat)
  #     }
  
  if (complete) {
    return(list(first = completeData, second = clinicalData))
    
  } else {
    return(list(first = unionData, second = clinicalData))
  }
}

#### Load in cases (full_data), not complete.
brca_full <- loadData(cancer = 'BRCA', clinicalDataBRCA, complete = TRUE)
kirc_full <- loadData(cancer = 'KIRC', clinicalDataKIRC, complete = TRUE)
lihc_full <- loadData(cancer = 'LIHC', clinicalDataLIHC, complete = TRUE)
luad_full <- loadData(cancer = 'LUAD', clinicalDataLUAD, complete = TRUE)

#### unlist kirc_full
kirc_data_full <- kirc_full[[1]]
brca_data_full <- brca_full[[1]]
lihc_data_full <- lihc_full[[1]]
luad_data_full <- luad_full[[1]]

kirc_clin <- kirc_full[[2]]
brca_clin <- brca_full[[2]]
lihc_clin <- lihc_full[[2]]
luad_clin <- luad_full[[2]]



# kirc_clin_full <- kirc_full[[2]]

####### transform ids
transform  <- function(data){
  transform_id_format <- function(x){
    x <- substr(x, 1, 12)
    x <- gsub('.', '-', x, fixed = TRUE)
    x <- tolower(x)
    
    return(x)
  }
  for(i in 1:3){
    colnames(data[[i]]) <- transform_id_format(colnames(data[[i]]))
  }
  return(data)
}

kirc_data_full <- transform(kirc_data_full)
brca_data_full <- transform(brca_data_full)
lihc_data_full <- transform(lihc_data_full)
luad_data_full <- transform(luad_data_full)


# Split in data types
kirc_methyl_full <- kirc_data_full[[1]]
kirc_mirna_full <- kirc_data_full[[2]]
kirc_mrna_full <- kirc_data_full[[3]]

brca_methyl_full <- brca_data_full[[1]]
brca_mirna_full <- brca_data_full[[2]]
brca_mrna_full <- brca_data_full[[3]]

lihc_methyl_full <- lihc_data_full[[1]]
lihc_mirna_full <- lihc_data_full[[2]]
lihc_mrna_full <- lihc_data_full[[3]]

luad_methyl_full <- luad_data_full[[1]]
luad_mirna_full <- luad_data_full[[2]]
luad_mrna_full <- luad_data_full[[3]]

#Merge data
dataMerge <- function(data, clinical){
  data <- as.data.frame(t(data))
  data <- cbind('id' = row.names(data), data)
  row.names(data) <- NULL
  names(clinical)[1] <- 'id'
  data <- inner_join(data, clinical, by = 'id')
  #data$days_to_death[is.na(as.numeric(data$days_to_death))] <- max(clinical$days_to_death, na.rm = T)
  return(data)
}

# factor_variables <- c(159, 161, 162, 163, 168, 170, 171, 176, 177, 198, 201, 206, 259, 260, 260, 261, 262, 
#                       271, 272, 278, 279, 283, 284, 285, 286, 288, 289, 295, 298, 309, 311, 312, 313, 315,
#                       316, 409, 418, 422, 426, 427, 428, 431, 432, 433, 436, 437, 438, 439)



pcaAll <- function(data, clin, name) {

  for (i in 7:8) {
    
    variable <- names(clin[i])
  
    
    ## choose clinical variable 
    # variable <- 'admin.project_code'
    
    temp <- clin[, c('bcr_patient_barcode', variable)]
    
    data_batch <- dataMerge(data,  temp)
    
   
    # kirc_mirna_merge_full <- dataMerge(kirc_mirna_full, temp)
    # kirc_mrna_merge_full <- dataMerge(kirc_mrna_full, temp)
    
    # # recode clinical barcode to id so can merge with kirc_methyl
    # names(clin)[2] <- 'id'
    # 
    # # inner join clin and kirc_methyl_merge_full by id
    # batch_data <- inner_join(kirc_methyl_merge_full, clin, by = 'id')
    # batch_data <- batch_data[, -c(20922:20926)]
    
    # run pca
    pca <- function(data){
    #data <- data[!is.na(data$days_to_death),]
    # data <- data[data$days_to_death > 0,]
    data_length <- (ncol(data)-1) 
    pca <- prcomp(data[,2:data_length])
    return(pca)
    }
    
    data_batch_pca <- pca(data_batch)
    
    
    # plot batch effect 
    pcaPlot <- function(pca, data,name){
    plot(pca$x[,1], 
         pca$x[,2],
         xlab = 'PCA 1',
         ylab = 'PCA 2',
         cex = 1,
         main = name,
         pch = 16,
         #        xlim= c(min, max),
         #        ylim = c(min, max),
         col = as.factor(data[,variable])
    )
    abline(v = c(0,0),
           h = c(0,0))
      
      legend('bottomright',
             legend = unique(as.factor(data[,variable][!is.na(data[,variable])])),
             col=1:length(as.factor(data[, variable])),
             pch=16,
             cex = 0.7)
    
     }
    
    pcaPlot(data_batch_pca, data_batch, name = name)
    print(i)

  }

}

pdf('/home/benbrew/Desktop/batches.pdf')
# run function
pcaAll(kirc_methyl_full, kirc_clin, 'kirc_methyl')
pcaAll(kirc_mirna_full, kirc_clin, 'kirc_mirna')
pcaAll(kirc_mrna_full, kirc_clin, 'kirc_mrna')

pcaAll(brca_methyl_full, brca_clin, 'brca_methyl')
pcaAll(brca_mirna_full, brca_clin, 'brca_mirna')
pcaAll(brca_mrna_full, brca_clin, 'brca_mrna')

pcaAll(lihc_methyl_full, lihc_clin, 'lihc_methyl')
pcaAll(lihc_mirna_full, lihc_clin, 'lihc_mirna')
pcaAll(lihc_mrna_full, lihc_clin, 'lihc_mrna')

pcaAll(luad_methyl_full, luad_clin, 'luad_methyl')
pcaAll(luad_mirna_full, luad_clin, 'luad_mirna')
pcaAll(luad_mrna_full, luad_clin, 'luad_mrna')


dev.off()

############### Look at to different vital statuses
vital1 <- names(clin[283])
vital2 <- names(clin[436])


############### Look closet at factor variable that showed pattern in pca
gender <- names(clin[285])
stage <- names(clin[418])
tissue <- names(clin[433])
year_diagnosis <- names(clin[439])



pcaFac <- function(variable) { 


  ## choose clinical variable 
  # variable <- 'admin.project_code'
  
  temp <- clin[, c('patient.bcr_patient_barcode', variable)]
  
  kirc_batch <- dataMerge(kirc_methyl_full,  temp)
  
  if (variable == year_diagnosis) { 
    kirc_batch$patient.year_of_initial_pathologic_diagnosis <- 
      ifelse(as.numeric(as.character(kirc_batch$patient.year_of_initial_pathologic_diagnosis)) > 2006, TRUE, FALSE)
  }
  
  # kirc_mirna_merge_full <- dataMerge(kirc_mirna_full, temp)
  # kirc_mrna_merge_full <- dataMerge(kirc_mrna_full, temp)
  
  # # recode clinical barcode to id so can merge with kirc_methyl
  # names(clin)[2] <- 'id'
  # 
  # # inner join clin and kirc_methyl_merge_full by id
  # batch_data <- inner_join(kirc_methyl_merge_full, clin, by = 'id')
  # batch_data <- batch_data[, -c(20922:20926)]
  
  # run pca
  pca <- function(data){
    #data <- data[!is.na(data$days_to_death),]
    # data <- data[data$days_to_death > 0,]
    data_length <- (ncol(data)-1) 
    pca <- prcomp(data[,2:data_length])
    return(pca)
  }
  
  batch_data_pca <- pca(kirc_batch)
  
  
  # plot batch effect 
  pcaPlot <- function(pca, data,name){
    plot(pca$x[,1], 
         pca$x[,2],
         xlab = 'PCA 1',
         ylab = 'PCA 2',
         cex = 1,
         main = name,
         pch = 16,
         #        xlim= c(min, max),
         #        ylim = c(min, max),
         col = as.factor(data[,variable])
    )
    abline(v = c(0,0),
           h = c(0,0))
    
    legend('bottomright',
           legend = unique(as.factor(data[,variable][!is.na(data[,variable])])),
           col=1:length(as.factor(data[, variable])),
           pch=16,
           cex = 0.7)
    
  }
  
  pcaPlot(batch_data_pca, kirc_batch, name = 'KIRC Methylation PCA')

}

pcaFac(gender)
pcaFac(year_diagnosis)
############### Numeric variables 

numeric_variables <- c(19, 20, 31, 32, 33, 35, 36, 156, 157, 158, 265, 267, 269, 270, 276, 297)


pdf('/home/benbrew/Desktop/batches_numeric.pdf')


for (i in numeric_variables) {
  
  variable <- names(clin[269])
  
  
  ## choose clinical variable 
  # variable <- 'admin.project_code'
  
  temp <- clin[, c('patient.bcr_patient_barcode', variable)]
  
  kirc_batch <- dataMerge(kirc_methyl_full,  temp)
  kirc_batch[, variable] <- as.numeric(kirc_batch[, variable])
  
  
  # kirc_mirna_merge_full <- dataMerge(kirc_mirna_full, temp)
  # kirc_mrna_merge_full <- dataMerge(kirc_mrna_full, temp)
  
  # # recode clinical barcode to id so can merge with kirc_methyl
  # names(clin)[2] <- 'id'
  # 
  # # inner join clin and kirc_methyl_merge_full by id
  # batch_data <- inner_join(kirc_methyl_merge_full, clin, by = 'id')
  # batch_data <- batch_data[, -c(20922:20926)]
  
  # run pca
  pca <- function(data){
    #data <- data[!is.na(data$days_to_death),]
    # data <- data[data$days_to_death > 0,]
    data_length <- (ncol(data)-1) 
    pca <- prcomp(data[,2:data_length])
    return(pca)
  }
  
  batch_data_pca <- pca(kirc_batch)
  
  
  # plot batch effect 
  pcaPlot <- function(pca, data,name){
    colVec<- colorRampPalette(c("green", "red"))(ceiling(max(data[, variable], na.rm = T)))
    data$cols <- colVec[ceiling(data[, variable])]
    
    plot(pca$x[,1], 
         pca$x[,2],
         xlab = 'PCA 1',
         ylab = 'PCA 2',
         cex = 1,
         main = name,
         pch = 16,
         #        xlim= c(min, max),
         #        ylim = c(min, max),
         col = data$cols
    )
    abline(v = c(0,0),
           h = c(0,0))
    
    
  }
  
  pcaPlot(batch_data_pca, kirc_batch, name = 'KIRC Methylation PCA')
  print(i)
  
}

dev.off()

#create a function for color vectors#####################################
pcaPlot <- function(pca, data, name){
  colVec<- colorRampPalette(c("green", "red"))(ceiling(max(data$days_to_death)))
  data$cols <- colVec[ceiling(data$days_to_death)]
  
  #   data$days_to_death <- (data$days_to_death - mean(data$days_to_death, na.rm = T))/
  #     sd(data$days_to_death, na.rm = T)
  min <- min(min(pca$x[,1]), pca$x[,2])
  max <- max(max(pca$x[,1]), pca$x[,2])
  plot(pca$x[,4], 
       pca$x[,5],
       xlab = 'PCA 1',
       ylab = 'PCA 2',
       cex = 1,
       main = name,
       pch = 16,
       xlim= c(min, max),
       ylim = c(min, max),
       col = adjustcolor(data$cols, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
}

################################################
# gender looks like it is the only interesting clinical variable




