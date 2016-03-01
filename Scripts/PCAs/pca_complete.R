#########
# This script will run a PCA on the data and plot it with patient health (days to death) as the size of the points

library(RColorBrewer)
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

######################################################################
# Load the original data
load_data <- function(cancer){

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

complete_ind <- feature_subset_indices(complete_data)
complete_data <- subset_data(complete_data, complete_ind)

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

complete_stat <- row_statistics(complete_data)
complete_data <- normalize_data(complete_data, complete_stat)
return(list(first = complete_data, second = clinical_data))

}

#### Load in cases txt files that have 2000 features and are normalized 
brca <- load_data(cancer = 'BRCA')
kirc <- load_data(cancer = 'KIRC')
lihc <- load_data(cancer = 'LIHC')
luad <- load_data(cancer = 'LUAD')


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
  colVec<- colorRampPalette(c("green", "red"))(ceiling(max(data$days_to_death)))
  data$cols <- colVec[ceiling(data$days_to_death)]
  
#   data$days_to_death <- (data$days_to_death - mean(data$days_to_death, na.rm = T))/
#     sd(data$days_to_death, na.rm = T)
  min <- min(min(pca$x[,1]), pca$x[,2])
  max <- max(max(pca$x[,1]), pca$x[,2])
  plot(pca$x[,4], 
       pca$x[,5],
       xlab = 'PCA 1',
       ylab = 'PCA 2',
       cex = 1,
       main = name,
       pch = 16,
       xlim= c(min, max),
       ylim = c(min, max),
       col = adjustcolor(data$cols, alpha.f = 0.5)
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

# library(ggbiplot)
# library("factoextra")
# 
# plot_vital <- function(pca, data, name){
#     fviz_pca_ind(pca,
#              geom = 'point', 
#              habillage = data$vital_status, #addEllipses=TRUE, ellipse.level=0.95,
#              cex = 5,
#              alpha = 0.6) + 
#       scale_color_brewer(palette="Set1") +
#       theme_minimal() +
#       labs(title = name)
# }
# 
# kirc_methyl$vital_status[1] <- 'alive'
# lihc_methyl$vital_status[1] <- 'alive'
# kirc_mirna$vital_status[1] <- 'alive'
# lihc_mirna$vital_status[1] <- 'alive'
# kirc_mrna$vital_status[1] <- 'alive'
# lihc_mrna$vital_status[1] <- 'alive'
# 
# par(mfrow = c(2,2))
# 
# plot_vital(brca_methyl_pca, brca_methyl, name = 'BRCA methyl')
# plot_vital(kirc_methyl_pca, kirc_methyl, name = 'KIRC methyl')
# plot_vital(lihc_methyl_pca, lihc_methyl, name = 'LIHC methyl')
# plot_vital(luad_methyl_pca, luad_methyl, name = 'LUAD methyl')
# 
# plot_vital(brca_mirna_pca, brca_mirna, 'BRCA mirna')
# plot_vital(kirc_mirna_pca, kirc_mirna, 'KIRC mirna')
# plot_vital(lihc_mirna_pca, lihc_mirna, 'LIHC mirna')
# plot_vital(luad_mirna_pca, luad_mirna, 'LUAD mirna')
# 
# plot_vital(brca_mrna_pca, brca_mrna, 'BRCA mrna')
# plot_vital(kirc_mrna_pca, kirc_mrna, 'KIRC mrna')
# plot_vital(lihc_mrna_pca, lihc_mrna, 'LIHC mrna')
# plot_vital(luad_mrna_pca, luad_mrna, 'LUAD mrna')
# 
# 
# # Load in labels for kirc and find outlier points 
# dataFolder <- paste(projectFolder, 'Scripts', '06_Two_Thousand_Features', 
#                     'cluster_complete_data',
#                     'Results', 'Labels', sep = '/')
# 
# # read in labels from complete data for hier, icluster and SNF
# hier <- read.table(paste0(dataFolder, '/2_1.txt'))
# icluster <- read.table(paste0(dataFolder, '/2_2.txt'))
# SNF <- read.table(paste0(dataFolder, '/2_3.txt'))
# 
# methods <- t(rbind(hier, icluster, SNF))
# 
# # Get complete IDs for each cancer type
# complete_ids <- colnames(kirc_methyl)
# kirc_labels <- as.data.frame(cbind(complete_ids, methods))
# row.names(kirc_labels) <- NULL
# colnames(kirc_labels) <- c('id', 'hier', 'icluster', 'SNF')
# 
# # Get position of outlier labels for pca methylation kircso
# temp <- sort(kirc_methyl_pca$x[,1])
# 
# # Positions are 171, 33, 38
# 
# 
# 
# 
