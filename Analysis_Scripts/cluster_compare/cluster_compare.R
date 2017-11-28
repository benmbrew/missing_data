# complete labels for cancers  
library(dplyr)
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "cluster_complete_data", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(testFolder, "Results/Labels", sep="/")

# hierarchical, icluster, snf 
# brca, kirc, lihc, luad

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

# load cluster complete data labels for each clustering method
kirc_hier <- as.data.frame(t(read.table(paste0(resultsFolder,'/2_1.txt'))))
kirc_iclus <- as.data.frame(t(read.table(paste0(resultsFolder,'/2_2.txt'))))
kirc_snf <- as.data.frame(t(read.table(paste0(resultsFolder,'/2_3.txt'))))

# add indicator to each label set 
kirc_hier$type <- 'hier'
kirc_iclus$type <- 'iclus'
kirc_snf$type <- 'snf'

# read in ids and bind them 
id_comp <- read.table('ids_complete')
kirc_hier <- cbind(kirc_hier, id_comp)
kirc_iclus <- cbind(kirc_iclus, id_comp)
kirc_snf <- cbind(kirc_snf, id_comp)

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
id_batch_data_450$id <- gsub('_', '-', id_batch_data_450$id)

kirc_hier[ ,3] <- transform_id_format(kirc_hier[, 3])
kirc_iclus[ ,3] <- transform_id_format(kirc_iclus[, 3])
kirc_snf[ ,3] <- transform_id_format(kirc_snf[, 3])

kirc_hier <- inner_join(kirc_hier, id_batch_data_450, by = c('x' = 'id'))
kirc_iclus <- inner_join(kirc_iclus, id_batch_data_450, by = c('x' = 'id'))
kirc_snf <- inner_join(kirc_snf, id_batch_data_450, by = c('x' = 'id'))


# group by gene and sum probe values. 
kirc_hier <- kirc_hier %>%
  group_by(V1) %>%
  summarise(batch_1275 = sum(batch == '1275'),
            batch_1418 = sum(batch == '1418'),
            batch_1424 = sum(batch == '1424'),
            batch_1500 = sum(batch == '1500'),
            batch_1536 = sum(batch == '1536'),
            batch_1670 = sum(batch == '1670'),
            batch_A264 = sum(batch == 'A264'),
            batch_A27A = sum(batch == 'A27A'),
            batch_A33L = sum(batch == 'A33L'),
            batch_A36Y= sum(batch == 'A36Y'),
            batch_A39G = sum(batch == 'A39G'))

# group by gene and sum probe values. 
kirc_snf <- kirc_snf %>%
  group_by(V1) %>%
  summarise(batch_1275 = sum(batch == '1275'),
            batch_1418 = sum(batch == '1418'),
            batch_1424 = sum(batch == '1424'),
            batch_1500 = sum(batch == '1500'),
            batch_1536 = sum(batch == '1536'),
            batch_1670 = sum(batch == '1670'),
            batch_A264 = sum(batch == 'A264'),
            batch_A27A = sum(batch == 'A27A'),
            batch_A33L = sum(batch == 'A33L'),
            batch_A36Y= sum(batch == 'A36Y'),
            batch_A39G = sum(batch == 'A39G'))



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



extractRelevantColumns <- function(data) {
  # List of features used for survival analysis
  features <- c("admin.batch_number",
                "patient.bcr_patient_barcode",
                "patient.bcr_patient_uuid",
                "patient.days_to_death",
                "patient.days_to_last_followup",
                "patient.days_to_last_known_alive",
                "patient.vital_status")
  #patientFeatures <- paste("patient", features, sep=".")
  
  # Add missing features to the data
  missingFeaturesInd <- !(features %in% colnames(data))
  data[features[missingFeaturesInd]] <- NA
  
  # Extract and rename the relevant columns
  data <- data[features]
  colnames(data) <- features
  
  return(data)
}

# Load clinical data with batch number
cancer <- 'KIRC'
# Process the clinical data
clin <- loadClinData(cancer)
clin <- extractRelevantColumns(clin)
