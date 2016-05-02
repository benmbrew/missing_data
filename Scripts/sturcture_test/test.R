################################################################################################
# This script tests if the missing data is random or structured
###############################################################################
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(testFolder, "Results", sep="/")

###############################################################################
# Initialize fixed variables
jvmGBLimit <- 8
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD","LUSC")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 12
numFeat <- 2000

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

# Load the original data
loadData <- function(cancer) {
  
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
  
  
  # transform patient IDs to the clinical ID format 
  transformIDFormat <- function(x){
    x <- substr(x, 1, 12)
    x <- gsub('.', '-', x, fixed = TRUE)
    x <- tolower(x)
    
    return(x)
  }
  
  
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
    clinicalIDs <- as.character(clinicalData$patient.bcr_patient_barcode)
    clinicalInd <- match(unionIDs, clinicalIDs)
    clinicalData <- clinicalData[clinicalInd, ]
  
    return(list(first = unionData, second = clinicalData))
}

#### Load in cases (full_data), not complete.
brca_full <- loadData(cancer = 'BRCA')
kirc_full <- loadData(cancer = 'KIRC')
lihc_full <- loadData(cancer = 'LIHC')
luad_full <- loadData(cancer = 'LUAD')
lusc_full <- loadData(cancer = 'LUSC')


####get gene data with three views
kirc_data_full <- kirc_full[[1]]
brca_data_full <- brca_full[[1]]
lihc_data_full <- lihc_full[[1]]
luad_data_full <- luad_full[[1]]
lusc_data_full <- lusc_full[[1]]

####get clinical data
kirc_clin <- kirc_full[[2]]
brca_clin <- brca_full[[2]]
lihc_clin <- lihc_full[[2]]
luad_clin <- luad_full[[2]]
lusc_clin <- lusc_full[[2]]


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
lusc_data_full <- transform(lusc_data_full)



# Split in data types
kirc_methyl <- kirc_data_full[[1]]
kirc_mirna <- kirc_data_full[[2]]
kirc_mrna <- kirc_data_full[[3]]

brca_methyl <- brca_data_full[[1]]
brca_mirna <- brca_data_full[[2]]
brca_mrna <- brca_data_full[[3]]

lihc_methyl <- lihc_data_full[[1]]
lihc_mirna <- lihc_data_full[[2]]
lihc_mrna <- lihc_data_full[[3]]

luad_methyl <- luad_data_full[[1]]
luad_mirna <- luad_data_full[[2]]
luad_mrna <- luad_data_full[[3]]

lusc_methyl <- lusc_data_full[[1]]
lusc_mirna <- lusc_data_full[[2]]
lusc_mrna <- lusc_data_full[[3]]

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

# 1) Could do t test or chi squared test
# 2) LittleMCAR on clinical data d 
# 3) logistic regression of each clincal variable on a 1/0 that indicates missingness.
# 4) MissMech (https://cran.r-project.org/web/packages/MissMech/MissMech.pdf)





