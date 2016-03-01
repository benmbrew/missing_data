############################################################
library(RColorBrewer)
library(impute)
library(dplyr)
setwd(dataFolder)

if('imputedDataPCA.RData' %in% dir()){
  load('imputedDataPCA.RData')
}else{
  
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "evaluate_original_imputation", sep="/")
resultsFolder <- paste(testFolder, "Results", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')

######################################################################
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
data_types <- c("methyl", "mirna", "mrna")

source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

# Load the original data
load_union_data <- function(cancer){
  
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
  
  # extract all cases which appear in all of the data types (intersection)
  union_data <- columnUnion(cases) # now ncol in each matrix is the same and identical to each other. 
  
  # subset the clinical data so that it corresponds to individuals in the union data
  union_ids <- colnames(union_data[[1]])
  union_ids <- transform_id_format(union_ids)
  
  # find the position of the patient IDS in the clinical data 
  clinical_data <- as.data.frame(clinical_data) # not in Daniel's original code 
  clinical_ids <- as.character(clinical_data$bcr_patient_barcode)
  clinical_ind <- match(union_ids, clinical_ids) # returns a vector of positions of (first) matches of its 
  # first argument in its second. Takes length of x with positions of y. NA where x is not in y. So this will
  # be length of union_ids with position of clinical ids where they match.
  clinical_data <- clinical_data[clinical_ind, ]# now clinical data has ids match with union data (cases)
  
  ######################################################################
  # Select a subset of features which differ most between cases and
  # controls.
  num_feat <- 2000
  num_views <- 3
  
  feature_subset_indices <- function(cases, subset_size = num_feat){
    num_views <- length(cases) # length of 3
    feature_subset_ind <- vector('list', num_views) # create vector with length of 3.
    
    for(v in 1:num_views){
      num_features <- nrow(cases[[v]])  
      pval <- sapply(1:num_features, 
                     function(i) t.test(cases[[v]][i,], 
                                        controls[[v]][i,])$p.value)  
      ind <- order(pval) 
      feature_subset_ind[[v]] <- ind[1:min(subset_size, num_features)] 
    }
    
    return(feature_subset_ind) 
  }
  
  subset_data <- function(data, ind){
    for(v in 1:length(data)){
      data[[v]] <- data[[v]][ind[[v]], ] # subsets data by the index that will keep the ones with the 
      # lowest pvals. These are the features that have the most significant difference between cases and controls.
      
    }
    return(data)
  }
  
  union_ind <- feature_subset_indices(union_data)
  union_data <- subset_data(union_data, union_ind)
  
  ######################################################################
  # Normalize the features in the data sets.
  # Normalization is performed before imputation and we expect that the
  # data will still be normalized after imputation (before clustering).
  row_statistics <- function(cases){
    
    num_views <- length(cases)
    row_stats <- vector('list', num_views)
    
    for(v in 1:num_views){
      #calculate the row means and std deviations 
      row_mean <- apply(cases[[v]], 1, mean, na.rm = T)
      row_sd <- apply(cases[[v]], 1, sd, na.rm = T)
      constant_ind <- row_sd == 0
      row_sd[constant_ind] <- 1
      row_stats[[v]] <- list(mean = row_mean, sd = row_sd, ind = constant_ind)
    }
    return(row_stats)
  }
  
  normalize_data <- function(data, stat){
    for(v in 1:length(data)) {
      data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
      data[[v]] <- data[[v]][!stat[[v]]$ind, ]
    }
    return(data)
  }
  
  union_stat <- row_statistics(union_data)
  union_data <- normalize_data(union_data, union_stat)
  return(list(first = union_data, second = clinical_data))
  
}

#### Load in cases txt files that have 2000 features and are normalized 
brca <- load_union_data(cancer = 'BRCA')
kirc <- load_union_data(cancer = 'KIRC')
lihc <- load_union_data(cancer = 'LIHC')
luad <- load_union_data(cancer = 'LUAD')

#### Get list of all 3 data types for each cancer. 
brca_data <- brca[[1]]

kirc_data <- kirc[[1]]

lihc_data <- lihc[[1]]

luad_data <- luad[[1]]

#### Use 4 imputation methods knnImputation, llsImputation, lsaImputation, randomImputation
imputeData <- function(data, sampleRows = FALSE){
  
  knn_data <- knnImputation(data, sampleRows)
  lls_data <- llsImputation(data, sampleRows)
  rand_data <- randomImputation(data, sampleRows)
  #lsa_kirc <- lsaImputation(kirc_data, sampleRows = FALSE)
  return(list(first = knn_data , second = lls_data, third = rand_data))
  
}

brca_impute <- imputeData(brca_data)
kirc_impute <- imputeData(kirc_data)
lihc_impute <- imputeData(lihc_data)
luad_impute <- imputeData(luad_data)

####### transform ids and unlist data
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


### BRCA
brca_knn <- brca_impute[[1]]
brca_knn <- transform(brca_knn)
brca_knn_methyl <- brca_knn[[1]]
brca_knn_mirna <- brca_knn[[2]]
brca_knn_mrna <- brca_knn[[3]]

brca_lls <- brca_impute[[2]]
brca_lls <- transform(brca_lls)
brca_lls_methyl <- brca_lls[[1]]
brca_lls_mirna <- brca_lls[[2]]
brca_lls_mrna <- brca_lls[[3]]

brca_rand <- brca_impute[[3]]
brca_rand <- transform(brca_rand)
brca_rand_methyl <- brca_rand[[1]]
brca_rand_mirna <- brca_rand[[2]]
brca_rand_mrna <- brca_rand[[3]]

### KIRC
kirc_knn <- kirc_impute[[1]]
kirc_knn <- transform(kirc_knn)
kirc_knn_methyl <- kirc_knn[[1]]
kirc_knn_mirna <- kirc_knn[[2]]
kirc_knn_mrna <- kirc_knn[[3]]

kirc_lls <- kirc_impute[[2]]
kirc_lls <- transform(kirc_lls)
kirc_lls_methyl <- kirc_lls[[1]]
kirc_lls_mirna <- kirc_lls[[2]]
kirc_lls_mrna <- kirc_lls[[3]]

kirc_rand <- kirc_impute[[3]]
kirc_rand <- transform(kirc_rand)
kirc_rand_methyl <- kirc_rand[[1]]
kirc_rand_mirna <- kirc_rand[[2]]
kirc_rand_mrna <- kirc_rand[[3]]

### LIHC
lihc_knn <- lihc_impute[[1]]
lihc_knn <- transform(lihc_knn)
lihc_knn_methyl <- lihc_knn[[1]]
lihc_knn_mirna <- lihc_knn[[2]]
lihc_knn_mrna <- lihc_knn[[3]]

lihc_lls <- lihc_impute[[2]]
lihc_lls <- transform(lihc_lls)
lihc_lls_methyl <- lihc_lls[[1]]
lihc_lls_mirna <- lihc_lls[[2]]
lihc_lls_mrna <- lihc_lls[[3]]

lihc_rand <- lihc_impute[[3]]
lihc_rand <- transform(lihc_rand)
lihc_rand_methyl <- lihc_rand[[1]]
lihc_rand_mirna <- lihc_rand[[2]]
lihc_rand_mrna <- lihc_rand[[3]]

### LUAD
luad_knn <- luad_impute[[1]]
luad_knn <- transform(luad_knn)
luad_knn_methyl <- luad_knn[[1]]
luad_knn_mirna <- luad_knn[[2]]
luad_knn_mrna <- luad_knn[[3]]

luad_lls <- luad_impute[[2]]
luad_lls <- transform(luad_lls)
luad_lls_methyl <- luad_lls[[1]]
luad_lls_mirna <- luad_lls[[2]]
luad_lls_mrna <- luad_lls[[3]]

luad_rand <- luad_impute[[3]]
luad_rand <- transform(luad_rand)
luad_rand_methyl <- luad_rand[[1]]
luad_rand_mirna <- luad_rand[[2]]
luad_rand_mrna <- luad_rand[[3]]


####### unlist clinical
brca_clin <- brca[[2]]
kirc_clin <- kirc[[2]]
lihc_clin <- lihc[[2]]
luad_clin <- luad[[2]]

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

### BRCA
brca_knn_methyl_merge <- dataMerge(brca_knn_methyl, brca_clin)
brca_knn_mirna_merge <- dataMerge(brca_knn_mirna, brca_clin)
brca_knn_mrna_merge <- dataMerge(brca_knn_mrna, brca_clin)

brca_lls_methyl_merge <- dataMerge(brca_lls_methyl, brca_clin)
brca_lls_mirna_merge <- dataMerge(brca_lls_mirna, brca_clin)
brca_lls_mrna_merge <- dataMerge(brca_lls_mrna, brca_clin)

brca_rand_methyl_merge <- dataMerge(brca_rand_methyl, brca_clin)
brca_rand_mirna_merge <- dataMerge(brca_rand_mirna, brca_clin)
brca_rand_mrna_merge <- dataMerge(brca_rand_mrna, brca_clin)

### KIRC
kirc_knn_methyl_merge <- dataMerge(kirc_knn_methyl, kirc_clin)
kirc_knn_mirna_merge <- dataMerge(kirc_knn_mirna, kirc_clin)
kirc_knn_mrna_merge <- dataMerge(kirc_knn_mrna, kirc_clin)

kirc_lls_methyl_merge <- dataMerge(kirc_lls_methyl, kirc_clin)
kirc_lls_mirna_merge <- dataMerge(kirc_lls_mirna, kirc_clin)
kirc_lls_mrna_merge <- dataMerge(kirc_lls_mrna, kirc_clin)

kirc_rand_methyl_merge <- dataMerge(kirc_rand_methyl, kirc_clin)
kirc_rand_mirna_merge <- dataMerge(kirc_rand_mirna, kirc_clin)
kirc_rand_mrna_merge <- dataMerge(kirc_rand_mrna, kirc_clin)

### LIHC
lihc_knn_methyl_merge <- dataMerge(lihc_knn_methyl, lihc_clin)
lihc_knn_mirna_merge <- dataMerge(lihc_knn_mirna, lihc_clin)
lihc_knn_mrna_merge <- dataMerge(lihc_knn_mrna, lihc_clin)

lihc_lls_methyl_merge <- dataMerge(lihc_lls_methyl, lihc_clin)
lihc_lls_mirna_merge <- dataMerge(lihc_lls_mirna, lihc_clin)
lihc_lls_mrna_merge <- dataMerge(lihc_lls_mrna, lihc_clin)

lihc_rand_methyl_merge <- dataMerge(lihc_rand_methyl, lihc_clin)
lihc_rand_mirna_merge <- dataMerge(lihc_rand_mirna, lihc_clin)
lihc_rand_mrna_merge <- dataMerge(lihc_rand_mrna, lihc_clin)

### LUAD
luad_knn_methyl_merge <- dataMerge(luad_knn_methyl, luad_clin)
luad_knn_mirna_merge <- dataMerge(luad_knn_mirna, luad_clin)
luad_knn_mrna_merge <- dataMerge(luad_knn_mrna, luad_clin)

luad_lls_methyl_merge <- dataMerge(luad_lls_methyl, luad_clin)
luad_lls_mirna_merge <- dataMerge(luad_lls_mirna, luad_clin)
luad_lls_mrna_merge <- dataMerge(luad_lls_mrna, luad_clin)

luad_rand_methyl_merge <- dataMerge(luad_rand_methyl, luad_clin)
luad_rand_mirna_merge <- dataMerge(luad_rand_mirna, luad_clin)
luad_rand_mrna_merge <- dataMerge(luad_rand_mrna, luad_clin)


############
# run pca and plot
# in plot make point size the days to death and point color vital status.

pca <- function(data){
  data <- data[!is.na(data$days_to_death),]
  data <- data[data$days_to_death > 0,]
  data_length <- (ncol(data)-5) 
  pca <- prcomp(data[,2:data_length])
  return(pca)
}

### BRCA
brca_knn_methyl_merge_pca <- pca(brca_knn_methyl_merge)
brca_knn_mirna_merge_pca <- pca(brca_knn_mirna_merge)
brca_knn_mrna_merge_pca <- pca(brca_knn_mrna_merge)

brca_lls_methyl_merge_pca <- pca(brca_lls_methyl_merge)
brca_lls_mirna_merge_pca <- pca(brca_lls_mirna_merge)
brca_lls_mrna_merge_pca <- pca(brca_lls_mrna_merge)

brca_rand_methyl_merge_pca <- pca(brca_rand_methyl_merge)
brca_rand_mirna_merge_pca <- pca(brca_rand_mirna_merge)
brca_rand_mrna_merge_pca <- pca(brca_rand_mrna_merge)

### KIRC
kirc_knn_methyl_merge_pca <- pca(kirc_knn_methyl_merge)
kirc_knn_mirna_merge_pca <- pca(kirc_knn_mirna_merge)
kirc_knn_mrna_merge_pca <- pca(kirc_knn_mrna_merge)

kirc_lls_methyl_merge_pca <- pca(kirc_lls_methyl_merge)
kirc_lls_mirna_merge_pca <- pca(kirc_lls_mirna_merge)
kirc_lls_mrna_merge_pca <- pca(kirc_lls_mrna_merge)

kirc_rand_methyl_merge_pca <- pca(kirc_rand_methyl_merge)
kirc_rand_mirna_merge_pca <- pca(kirc_rand_mirna_merge)
kirc_rand_mrna_merge_pca <- pca(kirc_rand_mrna_merge)

### LIHC
lihc_knn_methyl_merge_pca <- pca(lihc_knn_methyl_merge)
lihc_knn_mirna_merge_pca <- pca(lihc_knn_mirna_merge)
lihc_knn_mrna_merge_pca <- pca(lihc_knn_mrna_merge)

lihc_lls_methyl_merge_pca <- pca(lihc_lls_methyl_merge)
lihc_lls_mirna_merge_pca <- pca(lihc_lls_mirna_merge)
lihc_lls_mrna_merge_pca <- pca(lihc_lls_mrna_merge)

lihc_rand_methyl_merge_pca <- pca(lihc_rand_methyl_merge)
lihc_rand_mirna_merge_pca <- pca(lihc_rand_mirna_merge)
lihc_rand_mrna_merge_pca <- pca(lihc_rand_mrna_merge)

### LUAD
luad_knn_methyl_merge_pca <- pca(luad_knn_methyl_merge)
luad_knn_mirna_merge_pca <- pca(luad_knn_mirna_merge)
luad_knn_mrna_merge_pca <- pca(luad_knn_mrna_merge)

luad_lls_methyl_merge_pca <- pca(luad_lls_methyl_merge)
luad_lls_mirna_merge_pca <- pca(luad_lls_mirna_merge)
luad_lls_mrna_merge_pca <- pca(luad_lls_mrna_merge)

luad_rand_methyl_merge_pca <- pca(luad_rand_methyl_merge)
luad_rand_mirna_merge_pca <- pca(luad_rand_mirna_merge)
luad_rand_mrna_merge_pca <- pca(luad_rand_mrna_merge)

#### pcaPlot takes a pca object, data frame, and name and returns
#### a plot with PCA 1 and PCA2 with the color of the points corresponding to the 
#### samples days to death, with missing values replaced with max of days to death for that cancer. 
pcaPlot <- function(pca, data, name){
  colVec<- colorRampPalette(c("green", "red"))(ceiling(max(data$days_to_death, na.rm = T)))
  new_data <- data[data$days_to_death != 0,]
  new_data$cols <- colVec[new_data$days_to_death]
  
  #   data$days_to_death <- (data$days_to_death - mean(data$days_to_death, na.rm = T))/
  #     sd(data$days_to_death, na.rm = T)
#   min <- min(min(pca$x[,1]), pca$x[,2])
#   max <- max(max(pca$x[,1]), pca$x[,2])
  plot(pca$x[,1], 
       pca$x[,2],
       xlab = 'PCA 1',
       ylab = 'PCA 2',
       cex = 1,
       cex.main = 2,
       main = name,
       pch = 16,
       cex.axis = 2,
#        xlim= c(min, max),
#        ylim = c(min, max),
       col = adjustcolor(new_data$cols, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
}

########## Plot results

### BRCA
par(mfrow = c(1,3))
pcaPlot(brca_knn_methyl_merge_pca, name = 'brca knn methyl', data = brca_knn_methyl_merge)
pcaPlot(brca_knn_mirna_merge_pca, name = 'brca knn mirna', data = brca_knn_mirna_merge)
pcaPlot(brca_knn_mrna_merge_pca, name = 'brca knn mrna', data = brca_knn_mrna_merge)

pcaPlot(brca_lls_methyl_merge_pca, name = 'brca lls methyl', data = brca_lls_methyl_merge)
pcaPlot(brca_lls_mirna_merge_pca, name = 'brca lls mirna', data = brca_lls_mirna_merge)
pcaPlot(brca_lls_mrna_merge_pca, name = 'brca lls mrna', data = brca_lls_mrna_merge)

pcaPlot(brca_rand_methyl_merge_pca, name = 'brca rand methyl', data = brca_rand_methyl_merge)
pcaPlot(brca_rand_mirna_merge_pca, name = 'brca rand mirna', data = brca_rand_mirna_merge)
pcaPlot(brca_rand_mrna_merge_pca, name = 'brca rand mrna', data = brca_rand_mrna_merge)

### KIRC
pcaPlot(kirc_knn_methyl_merge_pca, name = 'kirc knn methyl', data = kirc_knn_methyl_merge)
pcaPlot(kirc_knn_mirna_merge_pca, name = 'kirc knn mirna', data = kirc_knn_mirna_merge)
pcaPlot(kirc_knn_mrna_merge_pca, name = 'kirc knn mrna', data = kirc_knn_mrna_merge)

pcaPlot(kirc_lls_methyl_merge_pca, name = 'kirc lls methyl', data = kirc_lls_methyl_merge)
pcaPlot(kirc_lls_mirna_merge_pca, name = 'kirc lls mirna', data = kirc_lls_mirna_merge)
pcaPlot(kirc_lls_mrna_merge_pca, name = 'kirc lls mrna', data = kirc_lls_mrna_merge)

pcaPlot(kirc_rand_methyl_merge_pca, name = 'kirc rand methyl', data = kirc_rand_methyl_merge)
pcaPlot(kirc_rand_mirna_merge_pca, name = 'kirc rand mirna', data = kirc_rand_mirna_merge)
pcaPlot(kirc_rand_mrna_merge_pca, name = 'kirc rand mrna', data = kirc_rand_mrna_merge)

### LIHC
pcaPlot(lihc_knn_methyl_merge_pca, name = 'lihc knn methyl', data = lihc_knn_methyl_merge)
pcaPlot(lihc_knn_mirna_merge_pca, name = 'lihc knn mirna', data = lihc_knn_mirna_merge)
pcaPlot(lihc_knn_mrna_merge_pca, name = 'lihc knn mrna', data = lihc_knn_mrna_merge)

pcaPlot(lihc_lls_methyl_merge_pca, name = 'lihc lls methyl', data = lihc_lls_methyl_merge)
pcaPlot(lihc_lls_mirna_merge_pca, name = 'lihc lls mirna', data = lihc_lls_mirna_merge)
pcaPlot(lihc_lls_mrna_merge_pca, name = 'lihc lls mrna', data = lihc_lls_mrna_merge)

pcaPlot(lihc_rand_methyl_merge_pca, name = 'lihc rand methyl', data = lihc_rand_methyl_merge)
pcaPlot(lihc_rand_mirna_merge_pca, name = 'lihc rand mirna', data = lihc_rand_mirna_merge)
pcaPlot(lihc_rand_mrna_merge_pca, name = 'lihc rand mrna', data = lihc_rand_mrna_merge)

### LUAD
pcaPlot(luad_knn_methyl_merge_pca, name = 'luad knn methyl', data = luad_knn_methyl_merge)
pcaPlot(luad_knn_mirna_merge_pca, name = 'luad knn mirna', data = luad_knn_mirna_merge)
pcaPlot(luad_knn_mrna_merge_pca, name = 'luad knn mrna', data = luad_knn_mrna_merge)

pcaPlot(luad_lls_methyl_merge_pca, name = 'luad lls methyl', data = luad_lls_methyl_merge)
pcaPlot(luad_lls_mirna_merge_pca, name = 'luad lls mirna', data = luad_lls_mirna_merge)
pcaPlot(luad_lls_mrna_merge_pca, name = 'luad lls mrna', data = luad_lls_mrna_merge)

pcaPlot(luad_rand_methyl_merge_pca, name = 'luad rand methyl', data = luad_rand_methyl_merge)
pcaPlot(luad_rand_mirna_merge_pca, name = 'luad rand mirna', data = luad_rand_mirna_merge)
pcaPlot(luad_rand_mrna_merge_pca, name = 'luad rand mrna', data = luad_rand_mrna_merge)

save.image('imputedDataPCA.RData')
}

