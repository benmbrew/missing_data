#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to extract the biological data for each patient

######################################################################
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben/"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
rawDataFolder <- paste(homeFolder, "projects/SNF/NM_2015/Data", sep="/")


######################################################################
# Initialize helper functions

loadData <- function(cancer, dataType, fileSuffix) {
  processingResult <- "BasicProcessingResult"
  rownamesIndex <- 1
  removeIndices <- 1:2
  
  if (dataType=="methyl") {
    processingResult <- "AdvancedProcessingResult"
  } else if (dataType=="mirna") {
    removeIndices <- 1
  } else if (dataType=="rna") {
    rownamesIndex <- 2
  }
  
  # Load the data
  data <- NULL
  fileName <- paste(cancer, dataType, fileSuffix, sep="_")
  filePath <- paste(rawDataFolder, cancer, processingResult, fileName,
                    sep="/")
  if (file.exists(filePath)) {
    data <- read.delim(filePath)
    
    # Rename the columns and rows
    colnames(data) <- substr(colnames(data), 1, 16)
    rownames(data) <- data[, rownamesIndex]
    
    # Convert the data frame to a numeric matrix
    data <- as.matrix(data[, -removeIndices])
  }
  
  return(data)
}

primarySolidTumors <- function(data) {
  tcgaIDs <- colnames(data)
  sampleType <- substr(tcgaIDs, 14, 15)
  return(sampleType == "01")
}

solidTissueNormals <- function(data) {
  tcgaIDs <- colnames(data)
  sampleType <- substr(tcgaIDs, 14, 15)
  return(sampleType == "11")
}

otherSampleTypes <- function(data) {
  tcgaIDs <- colnames(data)
  sampleType <- substr(tcgaIDs, 14, 15)
  return(sampleType[!(sampleType %in% c("01", "11"))])
}

missingFeatures <- function(data) {
  return(apply(data, 1, function(x) any(is.na(x))))
}

constantFeatures <- function(data) {
  return(apply(data, 1, function(x) var(x)==0))
}

logTransformData <- function(data) {
  return(log(data+1, base=2))
}

extractRelevantData <- function(data) {
  # Remove features which are missing from any sample
  data <- data[!missingFeatures(data), ]
  
  # Remove features with no variance
  data <- data[!constantFeatures(data), ]
  
  # Extract the tumor and normal samples
  cases <- data[, primarySolidTumors(data)]
  controls <- data[, solidTissueNormals(data)]
  
  return(list(cases=cases, controls=controls))
}

######################################################################
# Process and save the biological data

# Required Steps:
# -select the cancers to process
#   -BRCA, KIRC, LIHC, LUAD, LUSC
# -extract the three data types
#   -methyl, mirna, mrna
# -log transform the mirna and mrna data
# -remove missing features
# -remove constant features
# -separate the cases and controls
# -save the data using a concise names

cancers <- c("glioblastoma")
for (cancer in cancers) {
  # Process methylation
  methyl <- loadData(cancer, "methyl")
  methylList <- extractRelevantData(methyl)
  
  # Process miRNA
  # Need to check if both of these exist!
  mirna_ga <- loadData(cancer, "mirna", "ga__RPM.txt")
  mirna_hiseq <- loadData(cancer, "mirna", "hiseq__RPM.txt")
  mirna <- cbind(mirna_ga, mirna_hiseq)
  mirna <- logTransformData(mirna)
  mirnaList <- extractRelevantData(mirna)
  
  # Process gene expression
  mrna <- loadData(cancer, "rna", "gene.txt")
  mrna <- logTransformData(mrna)
  mrnaList <- extractRelevantData(mrna)
  
  casesList <- list(methylList[[1]], mirnaList[[1]], mrnaList[[1]])
  controlsList <- list(methylList[[2]], mirnaList[[2]], mrnaList[[2]])
  dataLists <- list(casesList, controlsList)

  # Save the data
  sampleTypes <- c("cases", "controls")
  dataTypes <- c("methyl", "mirna", "mrna")
  for (i in 1:length(sampleTypes)) {
    for (j in 1:length(dataTypes)) {
      fileName <- paste(cancer, "_", dataTypes[j], "_",
                        sampleTypes[i], ".txt", sep="")
      filePath <- paste(projectFolder, "Data", fileName, sep="/")
      write.table(dataLists[[i]][[j]], filePath, sep="\t")
    }
  }
}
