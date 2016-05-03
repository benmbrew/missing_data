################################################################################################
# This script tests if the missing data is random or structured
library(MASS)
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
    clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
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

# Add surv time to each clinical data 
survObject <- function(data) {
  
  data$survTime <- data$days_to_death
  
  # Replace missing survival times with days to last follow up
  missingSurvInd <- is.na(data$survTime)
  lastFollowup <- data$days_to_last_followup
  data$survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
  
  return(data)
  
}

brca_clin <- survObject(brca_clin)
kirc_clin <- survObject(kirc_clin)
lihc_clin <- survObject(lihc_clin)
luad_clin <- survObject(luad_clin)
lusc_clin <- survObject(lusc_clin)

###########################################################################################################
# Try to deal with structured missing data first with a t test or chi squared test 
# 1) Could do t test or chi squared test- create dummy variables for whether a variable is missing.
# 
# 1 = missing
# 0 = observed
# 
# You can then run t-tests and chi-square tests between this variable and other variables 
# in the data set to see if the missingness on this variable is related to the values of other variables.

# Write function that adds a column in clinical data for TRUE if missing in that data type 

addMissing <- function(methyl, mirna, mrna, clin) {
  
  clin$methyl_missing <- apply(methyl, 2, function(x) all(is.na(x)))
  clin$mirna_missing <- apply(mirna, 2, function(x) all(is.na(x)))
  clin$mrna_missing <- apply(mrna, 2, function(x) all(is.na(x)))
  
  return(clin)
}

brca_clin <- addMissing(brca_methyl, brca_mirna, brca_mrna, brca_clin)
kirc_clin <- addMissing(kirc_methyl, kirc_mirna, kirc_mrna, kirc_clin)
lihc_clin <- addMissing(lihc_methyl, lihc_mirna, lihc_mrna, lihc_clin)
luad_clin <- addMissing(luad_methyl, luad_mirna, luad_mrna, luad_clin)
lusc_clin <- addMissing(lusc_methyl, lusc_mirna, lusc_mrna, lusc_clin)

# Run t test on missing variables and each varible in clinical data (including survTime)

tTest <- function(data, data_type, column) {
  data <- data[rowSums(is.na(data)) < 7,]
  
  t.test(x = data[,column][which(data[,data_type] == "TRUE")],
         y = data[,column][which(data[,data_type] == "FALSE")])
}

# As the p-value 0.4828 is greater than the .05 significance level, we do not reject the null 
# hypothesis that the smoking habit is independent of the exercise level of the students.
cTest <- function(data, data_type, column) {
  
  data <- data[rowSums(is.na(data)) < 7,]
  

  tbl = table(data[,data_type], data[,column]) 
  chisq.test(tbl)
  
}
  
#### brca
tTest(brca_clin, "methyl_missing", "days_to_death")
tTest(brca_clin, "mirna_missing", "days_to_death")
tTest(brca_clin, "mrna_missing", "days_to_death")

tTest(brca_clin, "methyl_missing", "days_to_last_followup")
tTest(brca_clin, "mirna_missing", "days_to_last_followup")
tTest(brca_clin, "mrna_missing", "days_to_last_followup")

tTest(brca_clin, "methyl_missing", "survTime")
tTest(brca_clin, "mirna_missing", "survTime")
tTest(brca_clin, "mrna_missing", "survTime")

cTest(brca_clin, "methyl_missing", "vital_status")
cTest(brca_clin, "mirna_missing", "vital_status")
cTest(brca_clin, "mrna_missing", "vital_status")

#### kirc
tTest(kirc_clin, "methyl_missing", "days_to_death")
tTest(kirc_clin, "mirna_missing", "days_to_death")
tTest(kirc_clin, "mrna_missing", "days_to_death")

tTest(kirc_clin, "methyl_missing", "days_to_last_followup")
tTest(kirc_clin, "mirna_missing", "days_to_last_followup")
tTest(kirc_clin, "mrna_missing", "days_to_last_followup")

tTest(kirc_clin, "methyl_missing", "survTime")
tTest(kirc_clin, "mirna_missing", "survTime")
tTest(kirc_clin, "mrna_missing", "survTime")

cTest(kirc_clin, "methyl_missing", "vital_status")
cTest(kirc_clin, "mirna_missing", "vital_status")
cTest(kirc_clin, "mrna_missing", "vital_status")

#### lihc
tTest(lihc_clin, "methyl_missing", "days_to_death")
tTest(lihc_clin, "mirna_missing", "days_to_death")
tTest(lihc_clin, "mrna_missing", "days_to_death")

tTest(lihc_clin, "methyl_missing", "days_to_last_followup")
tTest(lihc_clin, "mirna_missing", "days_to_last_followup")
tTest(lihc_clin, "mrna_missing", "days_to_last_followup")

tTest(lihc_clin, "methyl_missing", "survTime")
tTest(lihc_clin, "mirna_missing", "survTime")
tTest(lihc_clin, "mrna_missing", "survTime")

cTest(lihc_clin, "methyl_missing", "vital_status")
cTest(lihc_clin, "mirna_missing", "vital_status")
cTest(lihc_clin, "mrna_missing", "vital_status")

#### luad
tTest(luad_clin, "methyl_missing", "days_to_death")
tTest(luad_clin, "mirna_missing", "days_to_death")
tTest(luad_clin, "mrna_missing", "days_to_death")

tTest(luad_clin, "methyl_missing", "days_to_last_followup")
tTest(luad_clin, "mirna_missing", "days_to_last_followup")
tTest(luad_clin, "mrna_missing", "days_to_last_followup")

tTest(luad_clin, "methyl_missing", "survTime")
tTest(luad_clin, "mirna_missing", "survTime")
tTest(luad_clin, "mrna_missing", "survTime")

cTest(luad_clin, "methyl_missing", "vital_status")
cTest(luad_clin, "mirna_missing", "vital_status")
cTest(luad_clin, "mrna_missing", "vital_status")

#### lusc
tTest(lusc_clin, "methyl_missing", "days_to_death")
tTest(lusc_clin, "mirna_missing", "days_to_death")
tTest(lusc_clin, "mrna_missing", "days_to_death")

tTest(lusc_clin, "methyl_missing", "days_to_last_followup")
tTest(lusc_clin, "mirna_missing", "days_to_last_followup")
tTest(lusc_clin, "mrna_missing", "days_to_last_followup")

tTest(lusc_clin, "methyl_missing", "survTime")
tTest(lusc_clin, "mirna_missing", "survTime")
tTest(lusc_clin, "mrna_missing", "survTime")

cTest(lusc_clin, "methyl_missing", "vital_status")
cTest(lusc_clin, "mirna_missing", "vital_status")
cTest(lusc_clin, "mrna_missing", "vital_status")

############################################################################################################
# 2) LittleMCAR on clinical data 


############################################################################################################
# 3) logistic regression of each clincal variable on a 1/0 that indicates missingness.

#############################################################################################################
# 4) MissMech (https://cran.r-project.org/web/packages/MissMech/MissMech.pdf)

