####################### Script for heat maps of the intersection
# One combined and 3 individual. 
library(SNFtool)
library(reshape2)
library(ggplot2)
library(gplots)
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
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

##################################################################

### function to take a list of data and create a spectral clustering graph of an SNF matrix
heatMap <- function(data_list, name){

  # change features to columns and samples to rows by transposing matrix.
  data <- lapply(data_list, t)
  
  # calculate distance between samples (each row)
  distances <- lapply(data, function(x) as.matrix(dist(x)))
  
  # convert the distances to  affinities
  affinities <- lapply(distances, affinityMatrix)
  
  # fuse the matrices 
  fused_matrix <- SNF(affinities)
  
  data_heatmap <- heatmap(fused_matrix, 
                          Rowv=NA, 
                          Colv=NA, 
                          labRow = FALSE,
                          labCol = FALSE,
                          col = heat.colors(256), 
                          scale="column", 
                          margins=c(5,10),
                          main = name)
  
#   # Choose number of clusters
#   num_clus_estimates <- unlist(estimateNumberOfClustersGivenGraph(fused_matrix, 2:10))
#   C <- max(num_clus_estimates[c(1,3)])
#   
#   group = spectralClustering(fused_matrix, C)
#   ## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.
#   
#   displayClusters(fused_matrix, group)
  
}

heatMap(brca_data, name = 'BRCA Intersection')
heatMap(kirc_data, name = 'KIRC Intersection')
heatMap(lihc_data, name = 'LIHC Intersection')
heatMap(luad_data, name = 'LUAD Intersection')


heatMapIndividual <- function(data, name, name2){ 
  data <- t(data)
  data_heatmap <- heatmap.2(data, dendrogram="none", 
                            Rowv=FALSE, 
                            Colv=FALSE, 
                            col = heat.colors(250), 
                            scale="none", 
                            key=TRUE, 
                            density.info="none", 
                            trace="none", 
                            cexRow=0.125, 
                            cexCol=0.125, 
                            symm=F,
                            symkey=T,
                            symbreaks=T,
                            main = name,
                            xlab = name2,
                            ylab = 'Samples')
                          
}

heatMapIndividual(brca_methyl, name = 'BRCA Methyl', name2 = 'methylation')
heatMapIndividual(kirc_methyl, name = 'KIRC Methyl', name2 = 'methylation')
heatMapIndividual(lihc_methyl, name = 'LIHC Methyl', name2 = 'methylation')
heatMapIndividual(luad_methyl, name = 'LUAD Methyl', name2 = 'methylation')

heatMapIndividual(brca_mirna, name = 'BRCA mirna', name2 = 'mirna')
heatMapIndividual(kirc_mirna, name = 'KIRC mirna', name2 = 'mirna')
heatMapIndividual(lihc_mirna, name = 'LIHC mirna', name2 = 'mirna')
heatMapIndividual(luad_mirna, name = 'LUAD mirna', name2 = 'mirna')

heatMapIndividual(brca_mrna, name = 'BRCA mrna', name2 = 'mrna')
heatMapIndividual(kirc_mrna, name = 'KIRC mrna', name2 = 'mrna')
heatMapIndividual(lihc_mrna, name = 'LIHC mrna', name2 = 'mrna')
heatMapIndividual(luad_mrna, name = 'LUAD mrna', name2 = 'mrna')






