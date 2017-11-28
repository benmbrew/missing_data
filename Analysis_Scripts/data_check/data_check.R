#############################
# PCA for intersection and union with all features. 
library(dplyr)

# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Scripts')
resultsFolder <- paste(testFolder, "Results", sep="/")

######################################################################
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
data_types <- c("methyl", "mirna", "mrna")

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

# Load the original data
load_full_data <- function(cancer, complete = FALSE){
  
  loadData <- function(dataType, suffix="") {
    fileName <- paste(cancer, "_", dataType, suffix,".txt", sep="")
    filePath <- paste(projectFolder, "Data", fileName, sep="/")
    return(read.delim(filePath))
  }
  
  num_views <- length(data_types)
  cases <- vector("list", num_views)
  controls <- vector("list", num_views)
  
  # Load the biological data
  for (v in 1:num_views) {
    cases[[v]] <- as.matrix(loadData(data_types[v], "_cases"))
    controls[[v]] <- as.matrix(loadData(data_types[v], "_controls"))
  }
  
  
  clinical_data <- loadData("clin")
  
  
  # transform patient IDs to the clinical ID format 
  transform_id_format <- function(x){
    x <- substr(x, 1, 12)
    x <- gsub('.', '-', x, fixed = TRUE)
    x <- tolower(x)
    
    return(x)
  }
  
  if(complete){
    # extract all cases which appear in all of the data types (intersection)
    complete_data <- columnIntersection(cases) # now ncol in each matrix is the same and identical to each other. 
    
    # subset the clinical data so that it corresponds to individuals in the complete data
    complete_ids <- colnames(complete_data[[1]])
    complete_ids <- transform_id_format(complete_ids)
    
    # find the position of the patient IDS in the clinical data 
    clinical_data <- as.data.frame(clinical_data) # not in Daniel's original code 
    clinical_ids <- as.character(clinical_data$bcr_patient_barcode)
    clinical_ind <- match(complete_ids, clinical_ids) # returns a vector of positions of (first) matches of its 
    # first argument in its second. Takes length of x with positions of y. NA where x is not in y. So this will
    # be length of complete_ids with position of clinical ids where they match.
    clinical_data <- clinical_data[clinical_ind, ]# now clinical data has ids match with complete data (cases)
  } else {
    
  }
  ######################################################################
  # Normalize the features in the data sets.
  # Normalization is performed before imputation and we expect that the
  #   # data will still be normalized after imputation (before clustering).
  #   row_statistics <- function(cases){
  #     
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
  #   normalize_data <- function(data, stat){
  #     for(v in 1:length(data)) {
  #       data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
  #       data[[v]] <- data[[v]][!stat[[v]]$ind, ]
  #     }
  #     return(data)
  #   }
  #   
  #   if(complete){
  #     completeStat <- row_statistics(complete_data)
  #     complete_data <- normalize_data(complete_data, completeStat)
  #   }else{
  #   incompleteStat <- row_statistics(cases)
  #   cases <- normalize_data(cases, incompleteStat)
  #   
  #   }
  if(complete){
    return(list(first = complete_data, second = clinical_data))
    # 
  }else{
    return(list(first = cases, second = clinical_data))
    
    
  }
  
}

#### Load in cases (full_data), not complete.
brca_full <- load_full_data(cancer = 'BRCA', complete = FALSE)
kirc_full <- load_full_data(cancer = 'KIRC', complete = FALSE)
lihc_full <- load_full_data(cancer = 'LIHC', complete = FALSE)
luad_full <- load_full_data(cancer = 'LUAD', complete = FALSE)

#### Load in complete data
brca_com <- load_full_data(cancer = 'BRCA', complete = TRUE)
kirc_com <- load_full_data(cancer = 'KIRC', complete = TRUE)
lihc_com <- load_full_data(cancer = 'LIHC', complete = TRUE)
luad_com <- load_full_data(cancer = 'LUAD', complete = TRUE)


### separate data types and clin data for full data
brca_data_full <- brca_full[[1]]
brca_clin_full <- brca_full[[2]]

kirc_data_full <- kirc_full[[1]]
kirc_clin_full <- kirc_full[[2]]

lihc_data_full <- lihc_full[[1]]
lihc_clin_full <- lihc_full[[2]]

luad_data_full <- luad_full[[1]]
luad_clin_full <- luad_full[[2]]

### separate data types and clin data for complete data
brca_data_com <- brca_com[[1]]
brca_clin_com <- brca_com[[2]]

kirc_data_com <- kirc_com[[1]]
kirc_clin_com <- kirc_com[[2]]

lihc_data_com <- lihc_com[[1]]
lihc_clin_com <- lihc_com[[2]]

luad_data_com <- luad_com[[1]]
luad_clin_com <- luad_com[[2]]

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

brca_data_full <- transform(brca_data_full)
kirc_data_full <- transform(kirc_data_full)
lihc_data_full <- transform(lihc_data_full)
luad_data_full <- transform(luad_data_full)

brca_data_com <- transform(brca_data_com)
kirc_data_com <- transform(kirc_data_com)
lihc_data_com <- transform(lihc_data_com)
luad_data_com <- transform(luad_data_com)


#### Split cases into methyl, mirna, and mrna
brca_methyl_full <- brca_data_full[[1]]
brca_mirna_full <- brca_data_full[[2]]
brca_mrna_full <- brca_data_full[[3]]

kirc_methyl_full <- kirc_data_full[[1]]
kirc_mirna_full <- kirc_data_full[[2]]
kirc_mrna_full <- kirc_data_full[[3]]

lihc_methyl_full <- lihc_data_full[[1]]
lihc_mirna_full <- lihc_data_full[[2]]
lihc_mrna_full <- lihc_data_full[[3]]

luad_methyl_full <- luad_data_full[[1]]
luad_mirna_full <- luad_data_full[[2]]
luad_mrna_full <- luad_data_full[[3]]

brca_methyl_com <- brca_data_com[[1]]
brca_mirna_com <- brca_data_com[[2]]
brca_mrna_com <- brca_data_com[[3]]

kirc_methyl_com <- kirc_data_com[[1]]
kirc_mirna_com <- kirc_data_com[[2]]
kirc_mrna_com <- kirc_data_com[[3]]

lihc_methyl_com <- lihc_data_com[[1]]
lihc_mirna_com <- lihc_data_com[[2]]
lihc_mrna_com <- lihc_data_com[[3]]

luad_methyl_com <- luad_data_com[[1]]
luad_mirna_com <- luad_data_com[[2]]
luad_mrna_com <- luad_data_com[[3]]

############################################
# Use column intersection and column union


############################################
# Look at complete data, number of '0.000000' in second and thrid view (mirna, mrna)
countZeros <- function(data) {
  return(length(which(data == 0))/(nrow(data)*ncol(data)))
}

countZeros(brca_mirna_com)
countZeros(brca_mrna_com)

countZeros(kirc_mirna_com)
countZeros(kirc_mrna_com)

countZeros(lihc_mirna_com)
countZeros(lihc_mrna_com)

countZeros(luad_mirna_com)
countZeros(luad_mrna_com)

# How many rows are entirely zero
rowZeroes <- function(data){
  count <- data[rowSums(data) == 0,]
  return(nrow(count)/nrow(data))
}
  
rowZeroes(brca_mirna_com)
rowZeroes(brca_mrna_com)

rowZeroes(kirc_mirna_com)
rowZeroes(kirc_mrna_com)

rowZeroes(lihc_mirna_com)
rowZeroes(lihc_mrna_com)

rowZeroes(luad_mirna_com)
rowZeroes(luad_mrna_com)


