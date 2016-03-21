
###############################################################################
# This script will look at the difference between clinical data in the intersection and union.

# Load libraries
library(dplyr)
library(ggplot2)
library(reshape2)

###############################################################################
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Combat",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(projectFolder, 'Scripts/06_Results', sep = '/')

###############################################################################
# Initialize fixed variables
jvmGBLimit <- 8
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 12
numFeat <- 2000

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))
source(paste0(resultsFolder, '/Lib/helpers.R'))

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
    
    
  } else {
  
    
    # Extract all cases which appear in at least one of the data types
    unionData <- columnUnion(cases)
    
    
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
#   if (complete){
#     completeInd <- featureSubsetIndices(completeData)
#     completeData <- subsetData(completeData, completeInd)
#   } else {
#     unionInd <- featureSubsetIndices(unionData)
#     unionData <- subsetData(unionData, unionInd)
#   }
#   #   
#   #   ####################################################################################
#   # #   Normalize the features in the data sets.
#   # #   Normalization is performed before imputation and we expect that the
# #   # #   data will still be normalized after imputation (before clustering).
#   rowStatistics <- function(cases){
#     num_views <- length(cases)
#     row_stats <- vector('list', num_views)
#     
#     for(v in 1:num_views){
#       #calculate the row means and std deviations 
#       row_mean <- apply(cases[[v]], 1, mean, na.rm = T)
#       row_sd <- apply(cases[[v]], 1, sd, na.rm = T)
#       constant_ind <- row_sd == 0
#       row_sd[constant_ind] <- 1
#       row_stats[[v]] <- list(mean = row_mean, sd = row_sd, ind = constant_ind)
#     }
#     return(row_stats)
#   }
#   
#   normalizeData <- function(data, stat){
#     for(v in 1:length(data)) {
#       data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
#       data[[v]] <- data[[v]][!stat[[v]]$ind, ]
#     }
#     return(data)
#   }
#   
#   if (complete) {
#     completeStat <- rowStatistics(completeData)
#     completeData <- normalizeData(completeData, completeStat)
#   } else {
#     unionStat <- rowStatistics(unionData)
#     unionData <- normalizeData(unionData, unionStat)
#   }
  
  if (complete) {
    return(completeData)
    
  } else {
    return(unionData)
  }
}

#####################################################################################
# Load complete, union and clinical data.

BRCAComplete <- loadData(cancer = 'BRCA', complete = TRUE)
BRCAUnion <- loadData(cancer = 'BRCA', complete = FALSE)

KIRCComplete <- loadData(cancer = 'KIRC', complete = TRUE)
KIRCUnion <- loadData(cancer = 'KIRC', complete = FALSE)

LIHCComplete <- loadData(cancer = 'LIHC', complete = TRUE)
LIHCUnion <- loadData(cancer = 'LIHC', complete = FALSE)

LUADComplete <- loadData(cancer = 'LUAD', complete = TRUE)
LUADUnion <- loadData(cancer = 'LUAD', complete = FALSE)

#######################################################################################
# separate into data types 
# BRCA
brca_comp_methyl <- BRCAComplete[[1]]
brca_comp_mirna <- BRCAComplete[[2]]
brca_comp_mrna <- BRCAComplete[[3]]

brca_union_methyl <- BRCAUnion[[1]]
brca_union_mirna <- BRCAUnion[[2]]
brca_union_mrna <- BRCAUnion[[3]]

# kirc
kirc_comp_methyl <- KIRCComplete[[1]]
kirc_comp_mirna <- KIRCComplete[[2]]
kirc_comp_mrna <- KIRCComplete[[3]]

kirc_union_methyl <- KIRCUnion[[1]]
kirc_union_mirna <- KIRCUnion[[2]]
kirc_union_mrna <- KIRCUnion[[3]]

# lihc
lihc_comp_methyl <- LIHCComplete[[1]]
lihc_comp_mirna <- LIHCComplete[[2]]
lihc_comp_mrna <- LIHCComplete[[3]]

lihc_union_methyl <- LIHCUnion[[1]]
lihc_union_mirna <- LIHCUnion[[2]]
lihc_union_mrna <- LIHCUnion[[3]]

# luad
luad_comp_methyl <- LUADComplete[[1]]
luad_comp_mirna <- LUADComplete[[2]]
luad_comp_mrna <- LUADComplete[[3]]

luad_union_methyl <- LUADUnion[[1]]
luad_union_mirna <- LUADUnion[[2]]
luad_union_mrna <- LUADUnion[[3]]

#######################################################################################
# extract samples in union that are not in intersection for each data type 

# BRCA
brca_not_methyl <- brca_union_methyl[,!colnames(brca_union_methyl) %in% colnames(brca_comp_methyl)]
brca_not_mirna <- brca_union_mirna[,!colnames(brca_union_mirna) %in% colnames(brca_comp_mirna)]
brca_not_mrna <- brca_union_mrna[,!colnames(brca_union_mrna) %in% colnames(brca_comp_mrna)]

# KIRC
kirc_not_methyl <- kirc_union_methyl[,!colnames(kirc_union_methyl) %in% colnames(kirc_comp_methyl)]
kirc_not_mirna <- kirc_union_mirna[,!colnames(kirc_union_mirna) %in% colnames(kirc_comp_mirna)]
kirc_not_mrna <- kirc_union_mrna[,!colnames(kirc_union_mrna) %in% colnames(kirc_comp_mrna)]

# lihc
lihc_not_methyl <- lihc_union_methyl[,!colnames(lihc_union_methyl) %in% colnames(lihc_comp_methyl)]
lihc_not_mirna <- lihc_union_mirna[,!colnames(lihc_union_mirna) %in% colnames(lihc_comp_mirna)]
lihc_not_mrna <- lihc_union_mrna[,!colnames(lihc_union_mrna) %in% colnames(lihc_comp_mrna)]

# luad
luad_not_methyl <- luad_union_methyl[,!colnames(luad_union_methyl) %in% colnames(luad_comp_methyl)]
luad_not_mirna <- luad_union_mirna[,!colnames(luad_union_mirna) %in% colnames(luad_comp_mirna)]
luad_not_mrna <- luad_union_mrna[,!colnames(luad_union_mrna) %in% colnames(luad_comp_mrna)]

###############################################################################################################
################# Compare kirc
# methyl
brca_comp_methyl_mean <- apply(brca_comp_methyl, 2, mean)
brca_not_methyl_mean <- apply(brca_not_methyl, 2, mean)
hist(brca_comp_methyl_mean, xlim = c(0, 1), col = adjustcolor('red', alpha.f = 0.6))
hist(brca_not_methyl_mean,xlim = c(0, 1), col = adjustcolor('blue', alpha.f = 0.6),add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))
# mirna
brca_comp_mirna_mean <- apply(brca_comp_mirna, 2, mean)
brca_not_mirna_mean <- apply(brca_not_mirna, 2, mean)
hist(brca_comp_mirna_mean, xlim = c(1, 3), col = adjustcolor('red', alpha.f = 0.6))
hist(brca_not_mirna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.6), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))


# mrna
brca_comp_mrna_mean <- apply(brca_comp_mrna, 2, mean)
brca_not_mrna_mean <- apply(brca_not_mrna, 2, mean)
hist(brca_not_mrna_mean, xlim = c(6, 7), col = adjustcolor('red', alpha.f = 0.5))
hist(brca_comp_mrna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.5), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))

################# Compare kirc
# methyl
kirc_comp_methyl_mean <- apply(kirc_comp_methyl, 2, mean)
kirc_not_methyl_mean <- apply(kirc_not_methyl, 2, mean)
hist(kirc_comp_methyl_mean, xlim = c(0.2, 0.6), col = adjustcolor('red', alpha.f = 0.6))
hist(kirc_not_methyl_mean,xlim = c(0, 1), col = adjustcolor('blue', alpha.f = 0.6),add = T)
legend('topleft', legend = c('complete', 'not complete'), fill = c('red', 'blue'))
# mirna
kirc_comp_mirna_mean <- apply(kirc_comp_mirna, 2, mean)
kirc_not_mirna_mean <- apply(kirc_not_mirna, 2, mean)
hist(kirc_comp_mirna_mean, xlim = c(1, 3), col = adjustcolor('red', alpha.f = 0.6))
hist(kirc_not_mirna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.6), add = T)
legend('topleft', legend = c('complete', 'not complete'), fill = c('red', 'blue'))


# mrna
kirc_comp_mrna_mean <- apply(kirc_comp_mrna, 2, mean)
kirc_not_mrna_mean <- apply(kirc_not_mrna, 2, mean)
hist(kirc_not_mrna_mean, xlim = c(5, 7.5),ylim = c(0, 150), col = adjustcolor('red', alpha.f = 0.5))
hist(kirc_comp_mrna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.5), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))

###############################################################################################################
################# Compare kirc
# methyl
brca_comp_methyl_mean <- apply(brca_comp_methyl, 2, mean)
brca_not_methyl_mean <- apply(brca_not_methyl, 2, mean)
hist(brca_comp_methyl_mean, xlim = c(0, 1), col = adjustcolor('red', alpha.f = 0.6))
hist(brca_not_methyl_mean,xlim = c(0, 1), col = adjustcolor('blue', alpha.f = 0.6),add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))
# mirna
brca_comp_mirna_mean <- apply(brca_comp_mirna, 2, mean)
brca_not_mirna_mean <- apply(brca_not_mirna, 2, mean)
hist(brca_comp_mirna_mean, xlim = c(1, 3), col = adjustcolor('red', alpha.f = 0.6))
hist(brca_not_mirna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.6), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))


# mrna
brca_comp_mrna_mean <- apply(brca_comp_mrna, 2, mean)
brca_not_mrna_mean <- apply(brca_not_mrna, 2, mean)
hist(brca_not_mrna_mean, xlim = c(6, 7), col = adjustcolor('red', alpha.f = 0.5))
hist(brca_comp_mrna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.5), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))

################# Compare lihc
# methyl
lihc_comp_methyl_mean <- apply(lihc_comp_methyl, 2, mean)
lihc_not_methyl_mean <- apply(lihc_not_methyl, 2, mean)
hist(lihc_comp_methyl_mean, xlim = c(0.2, 0.6), col = adjustcolor('red', alpha.f = 0.6))
hist(lihc_not_methyl_mean,xlim = c(0, 1), col = adjustcolor('blue', alpha.f = 0.6),add = T)
legend('topleft', legend = c('complete', 'not complete'), fill = c('red', 'blue'))
# mirna
lihc_comp_mirna_mean <- apply(lihc_comp_mirna, 2, mean)
lihc_not_mirna_mean <- apply(lihc_not_mirna, 2, mean)
hist(lihc_comp_mirna_mean, xlim = c(1, 3), col = adjustcolor('red', alpha.f = 0.6))
hist(lihc_not_mirna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.6), add = T)
legend('topleft', legend = c('complete', 'not complete'), fill = c('red', 'blue'))


# mrna
lihc_comp_mrna_mean <- apply(lihc_comp_mrna, 2, mean)
lihc_not_mrna_mean <- apply(lihc_not_mrna, 2, mean)
hist(lihc_not_mrna_mean, xlim = c(5, 7.5),ylim = c(0, 150), col = adjustcolor('red', alpha.f = 0.5))
hist(lihc_comp_mrna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.5), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))

###############################################################################################################
################# Compare kirc
# methyl
brca_comp_methyl_mean <- apply(brca_comp_methyl, 2, mean)
brca_not_methyl_mean <- apply(brca_not_methyl, 2, mean)
hist(brca_comp_methyl_mean, xlim = c(0, 1), col = adjustcolor('red', alpha.f = 0.6))
hist(brca_not_methyl_mean,xlim = c(0, 1), col = adjustcolor('blue', alpha.f = 0.6),add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))
# mirna
brca_comp_mirna_mean <- apply(brca_comp_mirna, 2, mean)
brca_not_mirna_mean <- apply(brca_not_mirna, 2, mean)
hist(brca_comp_mirna_mean, xlim = c(1, 3), col = adjustcolor('red', alpha.f = 0.6))
hist(brca_not_mirna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.6), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))


# mrna
brca_comp_mrna_mean <- apply(brca_comp_mrna, 2, mean)
brca_not_mrna_mean <- apply(brca_not_mrna, 2, mean)
hist(brca_not_mrna_mean, xlim = c(6, 7), col = adjustcolor('red', alpha.f = 0.5))
hist(brca_comp_mrna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.5), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))

################# Compare luad
# methyl
luad_comp_methyl_mean <- apply(luad_comp_methyl, 2, mean)
luad_not_methyl_mean <- apply(luad_not_methyl, 2, mean)
hist(luad_comp_methyl_mean, xlim = c(0.2, 0.6), col = adjustcolor('red', alpha.f = 0.6))
hist(luad_not_methyl_mean,xlim = c(0, 1), col = adjustcolor('blue', alpha.f = 0.6),add = T)
legend('topleft', legend = c('complete', 'not complete'), fill = c('red', 'blue'))
# mirna
luad_comp_mirna_mean <- apply(luad_comp_mirna, 2, mean)
luad_not_mirna_mean <- apply(luad_not_mirna, 2, mean)
hist(luad_comp_mirna_mean, xlim = c(1, 3), col = adjustcolor('red', alpha.f = 0.6))
hist(luad_not_mirna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.6), add = T)
legend('topleft', legend = c('complete', 'not complete'), fill = c('red', 'blue'))


# mrna
luad_comp_mrna_mean <- apply(luad_comp_mrna, 2, mean)
luad_not_mrna_mean <- apply(luad_not_mrna, 2, mean)
hist(luad_not_mrna_mean, xlim = c(5, 7.5),ylim = c(0, 150), col = adjustcolor('red', alpha.f = 0.5))
hist(luad_comp_mrna_mean, xlim = c(-3, 3), col = adjustcolor('blue', alpha.f = 0.5), add = T)
legend('topright', legend = c('complete', 'not complete'), fill = c('red', 'blue'))

#########################################################################################################
# look at genes for each cancer between complelte left out. 
#Methyl
brca_com_methyl_genes <- row.names(brca_comp_methyl)
brca_union_methyl_genes <- row.names(brca_union_methyl)
length(which(brca_com_methyl_genes %in% brca_union_methyl_genes))/length(brca_union_methyl_genes)
# 99%

kirc_com_methyl_genes <- row.names(kirc_comp_methyl)
kirc_union_methyl_genes <- row.names(kirc_union_methyl)
length(which(kirc_com_methyl_genes %in% kirc_union_methyl_genes))/length(kirc_union_methyl_genes)


lihc_com_methyl_genes <- row.names(lihc_comp_methyl)
lihc_union_methyl_genes <- row.names(lihc_union_methyl)
length(which(lihc_com_methyl_genes %in% lihc_union_methyl_genes))/length(lihc_union_methyl_genes)

luad_com_methyl_genes <- row.names(luad_comp_methyl)
luad_union_methyl_genes <- row.names(luad_union_methyl)
length(which(luad_com_methyl_genes %in% luad_union_methyl_genes))/length(luad_union_methyl_genes)


#mirna
brca_com_mirna_genes <- row.names(brca_comp_mirna)
brca_union_mirna_genes <- row.names(brca_union_mirna)
length(which(brca_com_mirna_genes %in% brca_union_mirna_genes))/length(brca_union_mirna_genes)


kirc_com_mirna_genes <- row.names(kirc_comp_mirna)
kirc_union_mirna_genes <- row.names(kirc_union_mirna)
length(which(kirc_com_mirna_genes %in% kirc_union_mirna_genes))/length(kirc_union_mirna_genes)


lihc_com_mirna_genes <- row.names(lihc_comp_mirna)
lihc_union_mirna_genes <- row.names(lihc_union_mirna)
length(which(lihc_com_mirna_genes %in% lihc_union_mirna_genes))/length(lihc_union_mirna_genes)


luad_com_mirna_genes <- row.names(luad_comp_mirna)
luad_union_mirna_genes <- row.names(luad_union_mirna)
length(which(luad_com_mirna_genes %in% luad_union_mirna_genes))/length(luad_union_mirna_genes)


#mrna
brca_com_mrna_genes <- row.names(brca_comp_mrna)
brca_union_mrna_genes <- row.names(brca_union_mrna)
length(which(brca_com_mrna_genes %in% brca_union_mrna_genes))/length(brca_union_mrna_genes)


kirc_com_mrna_genes <- row.names(kirc_comp_mrna)
kirc_union_mrna_genes <- row.names(kirc_union_mrna)
length(which(kirc_com_mrna_genes %in% kirc_union_mrna_genes))/length(kirc_union_mrna_genes)


lihc_com_mrna_genes <- row.names(lihc_comp_mrna)
lihc_union_mrna_genes <- row.names(lihc_union_mrna)
length(which(lihc_com_mrna_genes %in% lihc_union_mrna_genes))/length(lihc_union_mrna_genes)


luad_com_mrna_genes <- row.names(luad_comp_mrna)
luad_union_mrna_genes <- row.names(luad_union_mrna)
length(which(luad_com_mrna_genes %in% luad_union_mrna_genes))/length(luad_union_mrna_genes)


########################################################################################################
# Look at values of each gene in lihc for complete and union.

compareGene <- function(data, data2) {
  data <- t(data)
  data2 <- t(data2)
  data2 <- data2[complete.cases(data2),]
  data <- data[, colnames(data) %in% colnames(data2)]
  data2 <- data2[, colnames(data2) %in% colnames(data)]
  data_names <- colnames(data)
  
  scores <- matrix(, 0, 1)
  
  for(i in data_names)  {
    scores_temp <-  mean(data[,i], na.rm = T) - mean(data2[,i], na.rm = T)
    scores <- rbind(scores, scores_temp)
  }
  return(scores)
}

brca_methyl_compare <- compareGene(brca_comp_methyl, brca_not_methyl)
kirc_methyl_compare <- compareGene(kirc_comp_methyl, kirc_not_methyl)
lihc_methyl_compare <- compareGene(lihc_comp_methyl, lihc_not_methyl)
luad_methyl_compare <- compareGene(luad_comp_methyl, luad_not_methyl)
hist(brca_methyl_compare, ylim = c(0, 12000), xlim = c(-0.5,0.5), col = adjustcolor('red', alpha.f = 0.5))
hist(kirc_methyl_compare, ylim = c(0, 1100), xlim = c(0,2), col=adjustcolor('blue', alpha.f = 0.5), add = T)
hist(lihc_methyl_compare, ylim = c(0, 1100), xlim = c(0,2), col=adjustcolor('green', alpha.f = 0.5), add = T)
hist(luad_methyl_compare, ylim = c(0, 1100), xlim = c(0,2), col=adjustcolor('orange', alpha.f = 0.5), add = T)
legend('topright',legend = c('brca', 'kirc', 'lihc', 'luad'), fill = c('red', 'blue', 
                                                                      'green', 'orange'))


brca_mirna_compare <- compareGene(brca_comp_mirna, brca_not_mirna)
kirc_mirna_compare <- compareGene(kirc_comp_mirna, kirc_not_mirna)
lihc_mirna_compare <- compareGene(lihc_comp_mirna, lihc_not_mirna)
luad_mirna_compare <- compareGene(luad_comp_mirna, luad_not_mirna)
hist(brca_mirna_compare, ylim = c(0, 600), xlim = c(-1.5,1.5), col = adjustcolor('red', alpha.f = 0.5))
hist(kirc_mirna_compare, ylim = c(0, 600), xlim = c(0, 1.5), col=adjustcolor('blue', alpha.f = 0.5), add = T)
hist(lihc_mirna_compare, ylim = c(0, 600), xlim = c(0,1.5), col=adjustcolor('green', alpha.f = 0.5), add = T)
hist(luad_mirna_compare, ylim = c(0, 600), xlim = c(0,1.5), col=adjustcolor('orange', alpha.f = 0.5), add = T)
legend('topright',legend = c('brca', 'kirc', 'lihc', 'luad'), fill = c('red', 'blue', 
                                                                       'green', 'orange'))

brca_mrna_compare <- compareGene(brca_comp_mrna, brca_not_mrna)
kirc_mrna_compare <- compareGene(kirc_comp_mrna, kirc_not_mrna)
lihc_mrna_compare <- compareGene(lihc_comp_mrna, lihc_not_mrna)
luad_mrna_compare <- compareGene(luad_comp_mrna, luad_not_mrna)
hist(brca_mrna_compare, ylim = c(0, 15000), xlim = c(-1,1), col = adjustcolor('red', alpha.f = 0.5))
hist(kirc_mrna_compare, ylim = c(0, 1000), xlim = c(0, 0.5), col=adjustcolor('blue', alpha.f = 0.5), add = T)
hist(lihc_mrna_compare, ylim = c(0, 1000), xlim = c(0,0.5), col=adjustcolor('green', alpha.f = 0.5), add = T)
hist(luad_mrna_compare, ylim = c(0, 1000), xlim = c(0,0.5), col=adjustcolor('orange', alpha.f = 0.5), add = T)
legend('topright',legend = c('brca', 'kirc', 'lihc', 'luad'), fill = c('red', 'blue', 
                                                                       'green', 'orange'))






