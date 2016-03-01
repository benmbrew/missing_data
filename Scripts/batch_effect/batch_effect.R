###################################################
# This script is to examine the batch effect by looking at raw data and batch id
library(dplyr)
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(testFolder, "Results", sep="/")

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))


# Load methylation
# Load the data
cancer <- 'KIRC'
fileName <- paste(cancer, fileSuffix = 'methyl_450.txt', sep="_")
filePath <- paste(dataFolder, fileName,
                  sep="/")

methyl_450 <- read.delim(filePath, nrow = 10)


# extract the colum names from methyl_27
col_names_450 <- as.character(colnames(methyl_450))

################################################# for 450
#keep 1, 2, 3, 6 elements of each split
#split <- strsplit(sub_dat, 'X')
# last_split <- lapply(split, function(x) x[length(x)])
# colnames(data) <- unlist(last_split)
split_450 <- strsplit(col_names_450, '.', fixed = TRUE)
take_split_450 <- lapply(split_450, function(x) x[c(1:3,6 )])
data_frame_split_450 <- do.call('rbind', take_split_450)
id_batch_450 <- apply(data_frame_split_450[,1:3], 1, function(x) paste(x, collapse = '_'))
# now just split 1:3 together separate from 4
id_batch_data_450 <- cbind(id_batch_450, data_frame_split_450[, 4])
# id_batch_data_450[, 2] <- apply(id_batch_data_450[,2:4], 1, function(x) paste(x, collapse = '_'))
# id_batch_data_450 <- id_batch_data_450[, -c(3:4)]
id_batch_data_450 <- as.data.frame(id_batch_data_450)
names(id_batch_data_450) <- c('id', 'batch')
id_batch_data_450 <-id_batch_data_450[-1,]
rownames(id_batch_data_450) <- NULL

# drop duplicates and and make lowercase 
id_batch_data_450 <- id_batch_data_450[!duplicated(id_batch_data_450$id),]

# transform patient IDs to the clinical ID format 
transform_id_format <- function(x){
  x <- substr(x, 1, 12)
  x <- gsub('.', '-', x, fixed = TRUE)
  x <- tolower(x)
  
  return(x)
}

id_batch_data_450[,1] <- transform_id_format(id_batch_data_450[,1])



transposeDataFrame <- function(df, colnamesInd=1) {
  variableNames <- df[, colnamesInd]
  df <- as.data.frame(t(df[, -colnamesInd]))
  colnames(df) <- variableNames
  df
}

loadClinData <- function(cancer = 'KIRC') {
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


# 
# extractRelevantColumns <- function(data) {
#   # List of features used for survival analysis
#   features <- c("admin.batch_number",
#                 "patient.bcr_patient_barcode",
#                 "patient.bcr_patient_uuid",
#                 "patient.days_to_death",
#                 "patient.days_to_last_followup",
#                 "patient.days_to_last_known_alive",
#                 "patient.vital_status")
  #patientFeatures <- paste("patient", features, sep=".")
  
#   # Add missing features to the data
#   missingFeaturesInd <- !(features %in% colnames(data))
#   data[features[missingFeaturesInd]] <- NA
#   
#   # Extract and rename the relevant columns
#   data <- data[features]
#   colnames(data) <- features
#   
#   return(data)
# }

# Load clinical data with batch number
cancer <- 'KIRC'
# Process the clinical data
clin <- loadClinData(cancer)
# clin <- extractRelevantColumns(clin)

#################################
# Load kirc data 
data_types <- c("methyl", "mirna", "mrna")
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
kirc_full <- load_full_data(cancer = 'KIRC', complete = FALSE)

#### unlist kirc_full
kirc_data_full <- kirc_full[[1]]
kirc_clin_full <- kirc_full[[2]]

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

# Split in data types
kirc_methyl_full <- kirc_data_full[[1]]
kirc_mirna_full <- kirc_data_full[[2]]
kirc_mrna_full <- kirc_data_full[[3]]

#Merge data
dataMerge <- function(data, clinical){
  data <- as.data.frame(t(data))
  data <- cbind('id' = row.names(data), data)
  row.names(data) <- NULL
  names(clinical)[2] <- 'id'
  data <- inner_join(data, clinical, by = 'id')
  #data$days_to_death[is.na(as.numeric(data$days_to_death))] <- max(clinical$days_to_death, na.rm = T)
  return(data)
}

##### combine batch and clin 
id_batch_data_450$id <- gsub('_', '-', id_batch_data_450$id)

kirc_batch <- dataMerge(kirc_methyl_full,  clin)
kirc_mirna_merge_full <- dataMerge(kirc_mirna_full, temp)
kirc_mrna_merge_full <- dataMerge(kirc_mrna_full, temp)

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
  data_length <- (ncol(data)-6) 
  pca <- prcomp(data[,2:data_length])
  return(pca)
}

data_length <- (ncol(kirc_batch)-6) 
kirc_batch[, 2:data_length] <- scale(kirc_batch[, 2:data_length])

batch_data_pca <- pca(kirc_batch)
kirc_batch$col <- ifelse(grepl('A', kirc_batch$batch), 'green', 
                         ifelse(grepl('14', kirc_batch$batch), 'blue',
                                ifelse(grepl('15', kirc_batch$batch), 'black',
                                       ifelse(grepl('16', kirc_batch$batch), 'grey','yellow'))))

# plot batch effect 
pcaPlot <- function(pca, data,name){
  plot(pca$x[,3], 
       pca$x[,4],
       xlab = 'PCA 1',
       ylab = 'PCA 2',
       cex = 1,
       main = name,
       pch = 16,
       #        xlim= c(min, max),
       #        ylim = c(min, max),
       col = data$admin.batch_number
  )
  abline(v = c(0,0),
         h = c(0,0))
#   
#   legend('bottomright',
#          legend = unique(data$admin.batch_number),
#          col=1:length(data$admin.batch_number),
#          pch=16,
#          cex = 0.7)
  
}

pcaPlot(batch_data_pca, kirc_batch, name = 'batch')

temp <- kirc_batch[which(kirc_batch$admin.batch_number == '105.45.0' | kirc_batch$admin.batch_number == '63.53.0' | kirc_batch$admin.batch_number == '68.48.0' | 
                           kirc_batch$admin.batch_number == '70.45.0' | kirc_batch$admin.batch_number == '82.49.0' | kirc_batch$admin.batch_number == '90.41.0'),]

### look at characteristics of two clusters. pca$x[,2] < 0. 
batch_ind <-batch_data_pca$x[,2] < 0

kirc_batch <- kirc_batch[,20916:20921]

batch1 <- kirc_batch[batch_ind,]
batch2 <- kirc_batch[!batch_ind,]
