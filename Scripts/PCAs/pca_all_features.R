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

#########
#Merge data
dataMerge <- function(data, clinical){
  data <- as.data.frame(t(data))
  data <- cbind('id' = row.names(data), data)
  row.names(data) <- NULL
  names(clinical)[1] <- 'id'
  data <- left_join(data, clinical, by = 'id')
  data$days_to_death[is.na(as.numeric(data$days_to_death))] <- max(clinical$days_to_death, na.rm = T)
  return(data)
}

brca_methyl_merge_full <- dataMerge(brca_methyl_full, brca_clin_full)
brca_mirna_merge_full <- dataMerge(brca_mirna_full, brca_clin_full)
brca_mrna_merge_full <- dataMerge(brca_mrna_full, brca_clin_full)

kirc_methyl_merge_full <- dataMerge(kirc_methyl_full, kirc_clin_full)
kirc_mirna_merge_full <- dataMerge(kirc_mirna_full, kirc_clin_full)
kirc_mrna_merge_full <- dataMerge(kirc_mrna_full, kirc_clin_full)

lihc_methyl_merge_full <- dataMerge(lihc_methyl_full, lihc_clin_full)
lihc_mirna_merge_full <- dataMerge(lihc_mirna_full, lihc_clin_full)
lihc_mrna_merge_full <- dataMerge(lihc_mrna_full, lihc_clin_full)

luad_methyl_merge_full <- dataMerge(luad_methyl_full, luad_clin_full)
luad_mirna_merge_full <- dataMerge(luad_mirna_full, luad_clin_full)
luad_mrna_merge_full <- dataMerge(luad_mrna_full, luad_clin_full)

brca_methyl_merge_com <- dataMerge(brca_methyl_com, brca_clin_com)
brca_mirna_merge_com <- dataMerge(brca_mirna_com, brca_clin_com)
brca_mrna_merge_com <- dataMerge(brca_mrna_com, brca_clin_com)

kirc_methyl_merge_com <- dataMerge(kirc_methyl_com, kirc_clin_com)
kirc_mirna_merge_com <- dataMerge(kirc_mirna_com, kirc_clin_com)
kirc_mrna_merge_com <- dataMerge(kirc_mrna_com, kirc_clin_com)

lihc_methyl_merge_com <- dataMerge(lihc_methyl_com, lihc_clin_com)
lihc_mirna_merge_com <- dataMerge(lihc_mirna_com, lihc_clin_com)
lihc_mrna_merge_com <- dataMerge(lihc_mrna_com, lihc_clin_com)

luad_methyl_merge_com <- dataMerge(luad_methyl_com, luad_clin_com)
luad_mirna_merge_com <- dataMerge(luad_mirna_com, luad_clin_com)
luad_mrna_merge_com <- dataMerge(luad_mrna_com, luad_clin_com)

##### run pca on full data 
pca <- function(data){
  data <- data[!is.na(data$days_to_death),]
  data <- data[data$days_to_death > 0,]
  data_length <- (ncol(data)-5) 
  pca <- prcomp(data[,2:data_length])
  return(pca)
}

#### run pca on data
brca_methyl_full_pca <- pca(brca_methyl_merge_full)
brca_mirna_full_pca <- pca(brca_mirna_merge_full)
brca_mrna_full_pca <- pca(brca_mrna_merge_full)

kirc_methyl_full_pca <- pca(kirc_methyl_merge_full)
kirc_mirna_full_pca <- pca(kirc_mirna_merge_full)
kirc_mrna_full_pca <- pca(kirc_mrna_merge_full)

lihc_methyl_full_pca <- pca(lihc_methyl_merge_full)
lihc_mirna_full_pca <- pca(lihc_mirna_merge_full)
lihc_mrna_full_pca <- pca(lihc_mrna_merge_full)

luad_methyl_full_pca <- pca(luad_methyl_merge_full)
luad_mirna_full_pca <- pca(luad_mirna_merge_full)
luad_mrna_full_pca <- pca(luad_mrna_merge_full)


pcaPlot <- function(pca, data, data_com, name){
  data$color <- data$id %in% data_com$id
  data$color <- ifelse(data$color == TRUE, 'lightblue', 'black')
  #   data$days_to_death <- (data$days_to_death - mean(data$days_to_death, na.rm = T))/
  #     sd(data$days_to_death, na.rm = T)
  #   min <- min(min(pca$x[,1]), pca$x[,2])
  #   max <- max(max(pca$x[,1]), pca$x[,2])
  plot(pca$x[,1], 
       pca$x[,2],
       xlab = 'PCA 1',
       ylab = 'PCA 2',
       cex = 1,
       main = name,
       pch = 16,
       #        xlim= c(min, max),
       #        ylim = c(min, max),
       col = adjustcolor(data$color, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
}

########## Plot results
par(mfrow = c(2,2))
pcaPlot(brca_methyl_full_pca, name = 'BRCA methyl', brca_methyl_merge_full, brca_methyl_merge_com)
pcaPlot(kirc_methyl_full_pca, name = 'KIRC methyl', kirc_methyl_merge_full, kirc_methyl_merge_com)
pcaPlot(lihc_methyl_full_pca, name = 'LIHC methyl', lihc_methyl_merge_full, lihc_methyl_merge_com)
pcaPlot(luad_methyl_full_pca, name = 'LUAD methyl', luad_methyl_merge_full, luad_methyl_merge_com)

plot(brca_methyl_full_pca, type = 'l')
plot(kirc_methyl_full_pca, type = 'l')
plot(lihc_methyl_full_pca, type = 'l')
plot(luad_methyl_full_pca, type = 'l')

par(mfrow = c(2,2))
pcaPlot(brca_mirna_full_pca, name = 'BRCA mirna', brca_mirna_merge_full, brca_mirna_merge_com)
pcaPlot(kirc_mirna_full_pca, name = 'KIRC mirna', kirc_mirna_merge_full, kirc_mirna_merge_com)
pcaPlot(lihc_mirna_full_pca, name = 'LIHC mirna', lihc_mirna_merge_full, lihc_mirna_merge_com)
pcaPlot(luad_mirna_full_pca, name = 'LUAD mirna', luad_mirna_merge_full, luad_mirna_merge_com)

plot(brca_mirna_full_pca, type = 'l')
plot(kirc_mirna_full_pca, type = 'l')
plot(lihc_mirna_full_pca, type = 'l')
plot(luad_mirna_full_pca, type = 'l')


par(mfrow = c(2,2))
pcaPlot(brca_mrna_full_pca, name = 'BRCA mrna', brca_mrna_merge_full, brca_mrna_merge_com)
pcaPlot(kirc_mrna_full_pca, name = 'KIRC mrna', kirc_mrna_merge_full, kirc_mrna_merge_com)
pcaPlot(lihc_mrna_full_pca, name = 'LIHC mrna', lihc_mrna_merge_full, lihc_mrna_merge_com)
pcaPlot(luad_mrna_full_pca, name = 'LUAD mrna', luad_mrna_merge_full, luad_mrna_merge_com)

plot(brca_mrna_full_pca, type = 'l')
plot(kirc_mrna_full_pca, type = 'l')
plot(lihc_mrna_full_pca, type = 'l')
plot(luad_mrna_full_pca, type = 'l')

##########
# Look at kirc and get clinical data for the samples separated by the first PC. 
# kirc_methyl_merge_full$cluster_indicator <- kirc_methyl_full_pca$x[,1] > 0
# kirc_clinical <- kirc_methyl_merge_full[, 20916:20921]
# cluster1 <- kirc_clinical[kirc_clinical$cluster_indicator == TRUE,]
# cluster2 <- kirc_clinical[kirc_clinical$cluster_indicator == FALSE,]
# summary(cluster1)
# summary(cluster2)

#### Not much difference between two clusters in kirc methyl 

# ############
# ##### run pca on full data 
# 
# # Randomly remove 90% of data and see if clusters persist.
# data_length <- ncol(kirc_methyl_merge_full) - 6
# kirc_subset <- kirc_methyl_merge_full[, 2:data_length]
# kirc_subset <- kirc_subset[, sample(ncol(kirc_subset), ncol(kirc_subset)*.1) ]
# kirc_subset <- kirc_subset[]
# kirc_subset_pca <- prcomp(kirc_subset) 
# 
# plot(kirc_subset_pca$x[,1],
#            kirc_subset_pca$x[,2])
#####################################################################
# Look closer at kirc_methyl_merge

# column sums 
summary(colSums(kirc_methyl_full))

summary(colSums(brca_methyl_full))

summary(colSums(lihc_methyl_full))

summary(colSums(luad_methyl_full))

################################
summary(kirc_methyl_full_pca$x[,1])
summary(brca_methyl_full_pca$x[,1])
summary(lihc_methyl_full_pca$x[,1])
summary(luad_methyl_full_pca$x[,1])










