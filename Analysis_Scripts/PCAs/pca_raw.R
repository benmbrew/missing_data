##########################
# Load in raw data and do PCA of each cancer and data type without intersection or union.
library(dplyr)

# Initialize folders
homeFolder <- "/home/ben/Documents"
projectFolder <- paste0(homeFolder, "/missing_data")
testFolder <- paste0(projectFolder, "/Pipeline_Scripts")
dataFolder <- paste0(projectFolder, '/Data')
resultsFolder <- paste0(testFolder, "Results", sep="/")

######################################################################
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
data_types <- c("methyl", "mirna", "mrna")

source(paste(projectFolder, "Pipeline_Scripts/loadFunctions.R", sep="/"))

# Load the original data
load_raw_data <- function(cancer){
  
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
  
  incomplete_ind <- feature_subset_indices(cases)
  cases <- subset_data(cases, incomplete_ind)
  
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
  
  incompleteStat <- row_statistics(cases)
  cases <- normalize_data(cases, incompleteStat)
  return(list(first = cases, second = clinical_data))
  
}

#### Load in cases txt files that have 2000 features and are normalized 
brca <- load_raw_data(cancer = 'BRCA')
kirc <- load_raw_data(cancer = 'KIRC')
lihc <- load_raw_data(cancer = 'LIHC')
luad <- load_raw_data(cancer = 'LUAD')

### separate data types and clin data
brca_data <- brca[[1]]
brca_clin <- brca[[2]]

kirc_data <- kirc[[1]]
kirc_clin <- kirc[[2]]

lihc_data <- lihc[[1]]
lihc_clin <- lihc[[2]]

luad_data <- luad[[1]]
luad_clin <- luad[[2]]


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

brca_data <- transform(brca_data)
kirc_data <- transform(kirc_data)
lihc_data <- transform(lihc_data)
luad_data <- transform(luad_data)

#### Split cases into methyl, mirna, and mrna
brca_methyl <- brca_data[[1]]
brca_mirna <- brca_data[[2]]
brca_mrna <- brca_data[[3]]

kirc_methyl <- kirc_data[[1]]
kirc_mirna <- kirc_data[[2]]
kirc_mrna <- kirc_data[[3]]

lihc_methyl <- lihc_data[[1]]
lihc_mirna <- lihc_data[[2]]
lihc_mrna <- lihc_data[[3]]

luad_methyl <- luad_data[[1]]
luad_mirna <- luad_data[[2]]
luad_mrna <- luad_data[[3]]

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

brca_methyl_merge<- dataMerge(brca_methyl, brca_clin)
brca_mirna_merge<- dataMerge(brca_mirna, brca_clin)
brca_mrna_merge<- dataMerge(brca_mrna, brca_clin)

kirc_methyl_merge<- dataMerge(kirc_methyl, kirc_clin)
kirc_mirna_merge<- dataMerge(kirc_mirna, kirc_clin)
kirc_mrna_merge<- dataMerge(kirc_mrna, kirc_clin)

lihc_methyl_merge<- dataMerge(lihc_methyl, lihc_clin)
lihc_mirna_merge<- dataMerge(lihc_mirna, lihc_clin)
lihc_mrna_merge<- dataMerge(lihc_mrna, lihc_clin)

luad_methyl_merge<- dataMerge(luad_methyl, luad_clin)
luad_mirna_merge<- dataMerge(luad_mirna, luad_clin)
luad_mrna_merge<- dataMerge(luad_mrna, luad_clin)

###########
# use the intersected data from pca.R to find points in the intersection

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

#### run pca on data
brca_methyl_pca <- pca(brca_methyl_merge)
brca_mirna_pca <- pca(brca_mirna_merge)
brca_mrna_pca <- pca(brca_mrna_merge)

kirc_methyl_pca <- pca(kirc_methyl_merge)
kirc_mirna_pca <- pca(kirc_mirna_merge)
kirc_mrna_pca <- pca(kirc_mrna_merge)

lihc_methyl_pca <- pca(lihc_methyl_merge)
lihc_mirna_pca <- pca(lihc_mirna_merge)
lihc_mrna_pca <- pca(lihc_mrna_merge)

luad_methyl_pca <- pca(luad_methyl_merge)
luad_mirna_pca <- pca(luad_mirna_merge)
luad_mrna_pca <- pca(luad_mrna_merge)

#create a function for color vectors#####################################
pcaPlot <- function(pca, data, name){
  colVec<- colorRampPalette(c("green", "red"))(max(data$days_to_death))
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
       main = name,
       pch = 16,
#        xlim= c(min, max),
#        ylim = c(min, max),
       col = adjustcolor(new_data$cols, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
}

########## Plot results
par(mfrow = c(2,2))
pcaPlot(brca_methyl_pca, name = 'BRCA methyl', data = brca_methyl_merge)
pcaPlot(kirc_methyl_pca, name = 'KIRC methyl', data = kirc_methyl_merge)
pcaPlot(lihc_methyl_pca, name = 'LIHC methyl', data = lihc_methyl_merge)
pcaPlot(luad_methyl_pca, name = 'LUAD methyl', luad_methyl_merge)

plot(brca_methyl_pca, type = 'l')
plot(kirc_methyl_pca, type = 'l')
plot(lihc_methyl_pca, type = 'l')
plot(luad_methyl_pca, type = 'l')
# not as much variation accounted for in KIRC. 

pcaPlot(brca_mirna_pca, name = 'BRCA mirna', brca_mirna_merge)
pcaPlot(kirc_mirna_pca, name = 'KIRC mirna', kirc_mirna_merge)
pcaPlot(lihc_mirna_pca, name = 'LIHC mirna', lihc_mirna_merge)
pcaPlot(luad_mirna_pca, name = 'LUAD mirna', luad_mirna_merge)
plot(brca_mirna_pca, type = 'l')
plot(kirc_mirna_pca, type = 'l')
plot(lihc_mirna_pca, type = 'l')
plot(luad_mirna_pca, type = 'l')

pcaPlot(brca_mrna_pca, name  = 'BRCA mrna', brca_mrna_merge)
pcaPlot(kirc_mrna_pca, name  = 'KIRC mrna', kirc_mrna_merge)
pcaPlot(lihc_mrna_pca, name  = 'LIHC mrna', lihc_mrna_merge)
pcaPlot(luad_mrna_pca, name  = 'LUAD mrna', luad_mrna_merge)
plot(brca_mrna_pca, type = 'l')
plot(kirc_mrna_pca, type = 'l')
plot(lihc_mrna_pca, type = 'l')
plot(luad_mrna_pca, type = 'l')

### Plot with colored intersection
#create a function for color vectors#####################################
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
pcaPlot(brca_methyl_pca, name = 'BRCA methyl', brca_methyl_merge, brca_methyl_merge_com)
pcaPlot(kirc_methyl_pca, name = 'KIRC methyl', kirc_methyl_merge, kirc_methyl_merge_com)
pcaPlot(lihc_methyl_pca, name = 'LIHC methyl', lihc_methyl_merge, lihc_methyl_merge_com)
pcaPlot(luad_methyl_pca, name = 'LUAD methyl', luad_methyl_merge, luad_methyl_merge_com)

plot(brca_methyl_pca, type = 'l')
plot(kirc_methyl_pca, type = 'l')
plot(lihc_methyl_pca, type = 'l')
plot(luad_methyl_pca, type = 'l')
# not as much variation accounted for in KIRC. 
par(mfrow = c(2,2))
pcaPlot(brca_mirna_pca, name = 'BRCA mirna', brca_mirna_merge, brca_mirna_merge_com)
pcaPlot(kirc_mirna_pca, name = 'KIRC mirna', kirc_mirna_merge, kirc_mirna_merge_com)
pcaPlot(lihc_mirna_pca, name = 'LIHC mirna', lihc_mirna_merge, lihc_mirna_merge_com)
pcaPlot(luad_mirna_pca, name = 'LUAD mirna', luad_mirna_merge, luad_mirna_merge_com)

plot(brca_mirna_pca, type = 'l')
plot(kirc_mirna_pca, type = 'l')
plot(lihc_mirna_pca, type = 'l')
plot(luad_mirna_pca, type = 'l')


par(mfrow = c(2,2))
pcaPlot(brca_mrna_pca, name = 'BRCA mrna', brca_mrna_merge, brca_mrna_merge_com)
pcaPlot(kirc_mrna_pca, name = 'KIRC mrna', kirc_mrna_merge, kirc_mrna_merge_com)
pcaPlot(lihc_mrna_pca, name = 'LIHC mrna', lihc_mrna_merge, lihc_mrna_merge_com)
pcaPlot(luad_mrna_pca, name = 'LUAD mrna', luad_mrna_merge, luad_mrna_merge_com)

plot(brca_mrna_pca, type = 'l')
plot(kirc_mrna_pca, type = 'l')
plot(lihc_mrna_pca, type = 'l')
plot(luad_mrna_pca, type = 'l')
