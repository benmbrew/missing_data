################################################################################################
# This script tests if the missing data is random or structured
library(MASS)
library(BaylorEdPsych)
library(mvnmle)
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

### MAR
Run t test on missing variables and each varible in clinical data (including survTime)

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
# A large p-value (> 0.05) indicates weak evidence against the null hypothesis, 
# so you fail to reject the null hypothesis, in this case the null hypothesis is 
# that the data is MCAR, no patterns exists in the missing data.
# Little’s MCAR test is the most common test for missing cases being missing
# completely at random. If the p value for Little's MCAR test is not significant, then
# the data may be assumed to be MCAR and missingness is assumed not to matter
# for the analysis. 

# Missing completely at random (MCAR): data are missing independently of both
# observed and unobserved data. Example: a participant flips a coin to decide whether to complete the depression survey.
# Missing at random (MAR): given the observed data, data are missing
# independently of unobserved data. Example: male participants are more likely to refuse to fill out the depression
# survey, but it does not depend on the level of their depression.

# Roughly speaking, I think the "Little's test" you refer to looks to see if any observed variables are 
# predictive of the chance of other variables being missing, though I do not have the details at hand. 
# If there are variables predictive of non-response, the missingness mechanism is not MCAR but rather MAR, 
# and this needs to be kept in mind when doing the analysis (either by multiple imputation, EM algorithm or by 
# appropriate modelling) in order to ensure estimates and inferences are valid under MAR.

# Alternatively, Little (1988) provides a likelihood ratio test of
# the assumption of missing completely at random (MCAR). This test is part of
# the program BMDPAM, in the BMDP (Dixon, 1992) statistical package. In
# this example, the value of the likelihood ratio test is 2  256:61 p  0:000,
# indicating that the data in this example are not missing completely at random
# 
# A likelihood ratio test for MCAR exists. In statistics, a likelihood ratio test is a statistical 
# test used to compare the goodness of fit of two models, 
# one of which (the null model) is a special case of the other (the alternative model). 

# Add column for where TRUE is NA and false is 1 
LittleMcarMethyl <- function(data, vital_status = FALSE) {
  
  data$methyl_missing <- as.integer(ifelse(data$methyl_missing == TRUE, NA, 1))
  data$mirna_missing <-NULL
  data$mrna_missing <- NULL
  data <- data[rowSums(is.na(data)) < 8,]
  data <- data[!is.na(data$survTime),]
  data <- data[!is.na(data$vital_status),]
  data$vital_status <- data$vital_status
  missingInd <- is.na(data$methyl_missing)
  if (vital_status) {
    data$vital_status[missingInd] <- NA
    data <- data[, c(6,7)]
    
    LittleMCAR(data)
  } else {
    data$survTime[missingInd] <- NA
    data <- data[, c(6,7)]
    
    LittleMCAR(data)
    
  }
}

LittleMcarMirna <- function(data, vital_status = FALSE) {
  data$mirna_missing <- as.integer(ifelse(data$mirna_missing == TRUE, NA, 1))
  data$methyl_missing <-NULL
  data$mrna_missing <- NULL
  data <- data[rowSums(is.na(data)) < 8,]
  data <- data[!is.na(data$survTime),]
  data <- data[!is.na(data$vital_status),]
  data$vital_status <- data$vital_status
  missingInd <- is.na(data$mirna_missing)
  if (vital_status) {
    data$vital_status[missingInd] <- NA
    data <- data[, c(6,7)]
    
    LittleMCAR(data)
  } else {
    data$survTime[missingInd] <- NA
    data <- data[, c(6,7)]
    
    LittleMCAR(data)
    
  }
  
}

LittleMcarMrna <- function(data, vital_status = FALSE) {
  data$mrna_missing <- as.integer(ifelse(data$mrna_missing == TRUE, NA, 1))
  data$mirna_missing <-NULL
  data$methyl_missing <- NULL
  data <- data[rowSums(is.na(data)) < 8,]
  data <- data[!is.na(data$survTime),]
  data <- data[!is.na(data$vital_status),]
  data$vital_status <- data$vital_status
  missingInd <- is.na(data$mrna_missing)
  if (vital_status) {
    data$vital_status[missingInd] <- NA
    data <- data[, c(6,7)]
    
    LittleMCAR(data)
  } else {
    data$survTime[missingInd] <- NA
    data <- data[, c(6,7)]
    
    LittleMCAR(data)
    
  }
}

#BRCA
LittleMcarMethyl(brca_clin)$p.value
LittleMcarMethyl(brca_clin, vital_status = TRUE)$p.value
LittleMcarMirna(brca_clin)$p.value
LittleMcarMirna(brca_clin, vital_status = TRUE)$p.value
LittleMcarMrna(brca_clin)$p.value
LittleMcarMrna(brca_clin, vital_status = TRUE)$p.value

#KIRC
LittleMcarMethyl(kirc_clin)$p.value
LittleMcarMethyl(kirc_clin, vital_status = TRUE)$p.value
LittleMcarMirna(kirc_clin)$p.value
LittleMcarMirna(kirc_clin, vital_status = TRUE)$p.value
LittleMcarMrna(kirc_clin)$p.value
LittleMcarMrna(kirc_clin, vital_status = TRUE)$p.value

#LIHC
LittleMcarMethyl(lihc_clin)$p.value
LittleMcarMethyl(lihc_clin, vital_status = TRUE)$p.value
LittleMcarMirna(lihc_clin)$p.value
LittleMcarMirna(lihc_clin, vital_status = TRUE)$p.value
LittleMcarMrna(lihc_clin)$p.value
LittleMcarMrna(lihc_clin, vital_status = TRUE)$p.value

#LUAD
LittleMcarMethyl(luad_clin)$p.value
LittleMcarMethyl(luad_clin, vital_status = TRUE)$p.value
LittleMcarMirna(luad_clin)$p.value
LittleMcarMirna(luad_clin, vital_status = TRUE)$p.value
LittleMcarMrna(luad_clin)$p.value
LittleMcarMrna(luad_clin, vital_status = TRUE)$p.value

#LUSC
LittleMcarMethyl(lusc_clin)$p.value
LittleMcarMethyl(lusc_clin, vital_status = TRUE)$p.value
LittleMcarMirna(lusc_clin)$p.value
LittleMcarMirna(lusc_clin, vital_status = TRUE)$p.value
LittleMcarMrna(lusc_clin)$p.value
LittleMcarMrna(lusc_clin, vital_status = TRUE)$p.value


############################################################################################################
# 3) logistic regression of each clincal variable on a 1/0 that indicates missingness.
# The intercept is the predicted value of the dependent variable when all the independent variables are 0. 
# It is rarely of interest, unless the IVs are centered or standardized. Since you haven't 
# told us what any of your variables are, there's no way to say more than that.

#####BRCA
summary(glm(methyl_missing ~ days_to_death, data=brca_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_death, data=brca_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_death, data=brca_clin, family=binomial))
summary(glm(methyl_missing ~ days_to_last_followup, data=brca_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_last_followup, data=brca_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_last_followup, data=brca_clin, family=binomial))
summary(glm(methyl_missing ~ survTime, data=brca_clin, family=binomial))
summary(glm(mirna_missing ~ survTime, data=brca_clin, family=binomial))
summary(glm(mrna_missing ~ survTime, data=brca_clin, family=binomial))
summary(glm(methyl_missing ~ vital_status, data=brca_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status, data=brca_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status, data=brca_clin, family=binomial))

####kirc
summary(glm(methyl_missing ~ days_to_death, data=kirc_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_death, data=kirc_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_death, data=kirc_clin, family=binomial))
summary(glm(methyl_missing ~ days_to_last_followup, data=kirc_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_last_followup, data=kirc_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_last_followup, data=kirc_clin, family=binomial))
summary(glm(methyl_missing ~ survTime, data=kirc_clin, family=binomial))
summary(glm(mirna_missing ~ survTime, data=kirc_clin, family=binomial))
summary(glm(mrna_missing ~ survTime, data=kirc_clin, family=binomial))
summary(glm(methyl_missing ~ vital_status, data=kirc_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status, data=kirc_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status, data=kirc_clin, family=binomial))

####lihc
summary(glm(methyl_missing ~ days_to_death, data=lihc_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_death, data=lihc_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_death, data=lihc_clin, family=binomial))
summary(glm(methyl_missing ~ days_to_last_followup, data=lihc_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_last_followup, data=lihc_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_last_followup, data=lihc_clin, family=binomial))
summary(glm(methyl_missing ~ survTime, data=lihc_clin, family=binomial))
summary(glm(mirna_missing ~ survTime, data=lihc_clin, family=binomial))
summary(glm(mrna_missing ~ survTime, data=lihc_clin, family=binomial))
summary(glm(methyl_missing ~ vital_status, data=lihc_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status, data=lihc_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status, data=lihc_clin, family=binomial))

####luad
summary(glm(methyl_missing ~ days_to_death, data=luad_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_death, data=luad_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_death, data=luad_clin, family=binomial))
summary(glm(methyl_missing ~ days_to_last_followup, data=luad_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_last_followup, data=luad_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_last_followup, data=luad_clin, family=binomial))
summary(glm(methyl_missing ~ survTime, data=luad_clin, family=binomial))
summary(glm(mirna_missing ~ survTime, data=luad_clin, family=binomial))
summary(glm(mrna_missing ~ survTime, data=luad_clin, family=binomial))
summary(glm(methyl_missing ~ vital_status, data=luad_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status, data=luad_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status, data=luad_clin, family=binomial))


####lusc
summary(glm(methyl_missing ~ days_to_death, data=lusc_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_death, data=lusc_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_death, data=lusc_clin, family=binomial))
summary(glm(methyl_missing ~ days_to_last_followup, data=lusc_clin, family=binomial))
summary(glm(mirna_missing ~ days_to_last_followup, data=lusc_clin, family=binomial))
summary(glm(mrna_missing ~ days_to_last_followup, data=lusc_clin, family=binomial))
summary(glm(methyl_missing ~ survTime, data=lusc_clin, family=binomial))
summary(glm(mirna_missing ~ survTime, data=lusc_clin, family=binomial))
summary(glm(mrna_missing ~ survTime, data=lusc_clin, family=binomial))
summary(glm(methyl_missing ~ vital_status, data=lusc_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status, data=lusc_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status, data=lusc_clin, family=binomial))

############################################################################################################
# 3) Multiple logistic regression of each clincal variable on a 1/0 that indicates missingness.

#####BRCA
summary(glm(methyl_missing ~ vital_status + survTime, data=brca_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status + survTime, data=brca_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status + survTime, data=brca_clin, family=binomial))

####kirc
summary(glm(methyl_missing ~ vital_status + survTime, data=kirc_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status + survTime, data=kirc_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status + survTime, data=kirc_clin, family=binomial))

####lihc
summary(glm(methyl_missing ~ vital_status + survTime, data=lihc_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status + survTime, data=lihc_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status + survTime, data=lihc_clin, family=binomial))

####luad
summary(glm(methyl_missing ~ vital_status + survTime, data=luad_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status + survTime, data=luad_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status + survTime, data=luad_clin, family=binomial))


####lusc
summary(glm(methyl_missing ~ vital_status + survTime, data=lusc_clin, family=binomial))
summary(glm(mirna_missing ~ vital_status + survTime, data=lusc_clin, family=binomial))
summary(glm(mrna_missing ~ vital_status + survTime, data=lusc_clin, family=binomial))

