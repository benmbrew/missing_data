#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to extract the clinical data for each patient

######################################################################
# Initialize folders
homeFolder <- "/hpf/largeprojects/agoldenb/daniel"
rawDataFolder <- paste(homeFolder, "Data/TCGA/2014_12_06", sep="/")
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")

######################################################################
# Initialize helper functions

transposeDataFrame <- function(df, colnamesInd=1) {
  variableNames <- df[, colnamesInd]
  df <- as.data.frame(t(df[, -colnamesInd]))
  colnames(df) <- variableNames
  df
}

loadClinData <- function(cancer) {
  processingResult <- "RawData"
  fileSuffix <- "clinical.txt"
  
  # Load the data
  data <- NULL
  fileName <- paste(cancer, fileSuffix, sep="_")
  filePath <- paste(rawDataFolder, cancer, processingResult, fileName,
                    sep="/")
  if (file.exists(filePath)) {
    data <- read.delim(filePath)
    
    # Transpose the data
    data <- transposeDataFrame(data)
  }
  
  return(data)
}

extractRelevantColumns <- function(data) {
  # List of features used for survival analysis
  features <- c("bcr_patient_barcode",
                "bcr_patient_uuid",
                "days_to_death",
                "days_to_last_followup",
                "days_to_last_known_alive",
                "vital_status")
  patientFeatures <- paste("patient", features, sep=".")
  
  # Add missing features to the data
  missingFeaturesInd <- !(patientFeatures %in% colnames(data))
  data[patientFeatures[missingFeaturesInd]] <- NA
  
  # Extract and rename the relevant columns
  data <- data[patientFeatures]
  colnames(data) <- features
  
  return(data)
}

######################################################################
# Process and save the clinical data

cancers <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
for (cancer in cancers) {
  # Process the clinical data
  clin <- loadClinData(cancer)
  clin <- extractRelevantColumns(clin)
  
  # Save the data
  fileName <- paste(cancer, "clin.txt", sep="_")
  filePath <- paste(projectFolder, "Data", fileName, sep="/")
  write.table(clin, filePath, quote=FALSE, sep="\t", row.names=FALSE)
}