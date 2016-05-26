#########################################################################################################
# This script will compare labels from the intersection with those from the union. 

# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)
library(dplyr)
library(sva)

#########################################################################################################
# initiate folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
completeFolder <- paste(projectFolder, "Scripts",
                        "06_Cluster_Sizes",
                        "cluster_complete_data", 'Results', 'Labels',sep="/")
imputeOrigFolder <- paste(projectFolder, "Scripts",
                          "06_Cluster_Sizes",
                          "evaluate_original_imputation/Results/Clustering", sep="/")
similarityOrigFolder <- paste(projectFolder, "Scripts",
                              "06_Cluster_Sizes",
                              "evaluate_original_similarity/Results/Similarity", sep="/")
imputeFolder <- paste(projectFolder, "Scripts",
                      "06_Cluster_Sizes",
                      "evaluate_imputation/Results/Clustering", sep="/")
similarityFolder <- paste(projectFolder, "Scripts",
                          "06_Cluster_Sizes",
                          "evaluate_similarity/Results/Similarity", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Cluster_Compare", sep="/")
idsFolder <- paste(testFolder, "ids", sep="/")

#########################################################################################################
# Get ids from data for each cancer
source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))
dataTypes <- c("methyl", "mirna", "mrna")

# Load the original data
loadIDs <- function(cancer, combat = FALSE){
  
  if (cancer == "KIRC" & combat) {
    #####################################################################
    # Load in raw clinical data for kirc
    transposeDataFrame <- function(df, colnamesInd=1) {
      variableNames <- df[, colnamesInd]
      df <- as.data.frame(t(df[, -colnamesInd]))
      colnames(df) <- variableNames
      df
    }
    
    
    extractRelevantColumns <- function(data) {
      # List of features used for survival analysis
      features <- c("admin.batch_number",
                    "patient.bcr_patient_barcode",
                    "patient.bcr_patient_uuid",
                    "patient.days_to_death",
                    "patient.days_to_last_followup",
                    "patient.days_to_last_known_alive",
                    "patient.vital_status",
                    "patient.gender")
      patientFeatures <- paste("patient", features, sep=".")
      
      # Add missing features to the data
      missingFeaturesInd <- !(features %in% colnames(data))
      data[features[missingFeaturesInd]] <- NA
      
      # Extract and rename the relevant columns
      data <- data[features]
      colnames(data) <- features
      
      return(data)
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
    
    # Load clinical data with batch number
    cancer <- 'KIRC'
    # Process the clinical data
    clinicalData <- loadClinData(cancer)
    clinicalData <- extractRelevantColumns(clinicalData)
    
    ## remove patient from column names
    features <- colnames(clinicalData)
    split <- strsplit(features, '.', fixed = TRUE)
    keepSplit <- lapply(split, function(x) x[length(x)])
    features <- unlist(keepSplit)
    colnames(clinicalData) <- features
    
    ######################################################################
    # Load functions
    
    # Note: some functions depend on variables initialized above!
    # lsaImputation:
    # -imputedFile, incompleteFile, projectFolder, jvmGBLimit
    # iClusterClustering:
    # -numCores
    source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))
    
    ######################################################################
    # Load the original data
    
    loadData <- function(dataType, suffix="") {
      fileName <- paste(cancer, "_", dataType, suffix,".txt", sep="")
      filePath <- paste(projectFolder, "Data", fileName, sep="/")
      return(as.matrix(read.delim(filePath)))
    }
    
    numViews <- length(dataTypes)
    cases <- vector("list", numViews)
    controls <- vector("list", numViews)
    
    # Load the biological data
    for (v in 1:numViews) {
      cases[[v]] <- loadData(dataTypes[v], "_cases")
      controls[[v]] <- loadData(dataTypes[v], "_controls")
    }
    
    #####################################################################
    # Generate subsets of the cases derived from the set of individuals
    # which are present in all data types
    
    # Extract all cases which appear in all of the data types
    completeData <- columnIntersection(cases)
    
    #####################################################################
    # subset clinical data by ids in cases
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
    
    transformUnionId <- function(x) {
      # Keep the first 12 characters
      x <- substr(x, 1, 12)
      
      return(x)
    }
    
    unionData <- columnUnion(cases)
    
    # Subset the clinical data so that it corresponds to individuals
    # in the union data
    
    for (i in 1:3) {
      temp.data  <- unionData[[i]]
      temp.names <- transformUnionId(colnames(temp.data))
      colnames(temp.data) <- temp.names
      temp.data <- temp.data[, !duplicated(colnames(temp.data))]
      unionData[[i]] <- temp.data
    }
    
    # Subset the clinical data so that it corresponds to individuals
    # in the union data
    unionIDs <- colnames(unionData[[1]])
    unionIDs <- transformIDFormat(unionIDs)
    # Find the position of the patient IDs in the clinical data
    clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
    clinicalInd <- match(unionIDs, clinicalIDs)
    clinicalData <- clinicalData[clinicalInd, ]
    
    # Subset clinical data and completeData by clinicalData
    clinicalData <- clinicalData[rowSums(is.na(clinicalData)) < ncol(clinicalData),]
    clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
    unionInd <- match(clinicalIDs, unionIDs)
    
    clinicalData$days_to_death <- as.numeric(clinicalData$days_to_death)
    clinicalData$days_to_last_followup <- as.numeric(clinicalData$days_to_last_followup)
    ###### Run combat function 
    for (i in 1:numViews) {
      temp.unionData <- unionData[[i]]
      temp.unionData <- temp.unionData[,unionInd]
      temp.2.unionData <- temp.unionData[rowSums(temp.unionData, na.rm = T) != 0,]
      temp.modcombat <- model.matrix(~1, data = clinicalData)
      temp.batch <- clinicalData$gender
      temp_combat = ComBat(dat=temp.2.unionData, batch=temp.batch, mod=temp.modcombat, par.prior=TRUE, prior.plots=FALSE)
      unionData[[i]] <- temp_combat
    }
    
    unionIds <- colnames(unionData[[1]])
    
    
    completeIDs <- colnames(completeData[[1]])
    completeIDs <- transformIDFormat(completeIDs)
    # Find the position of the patient IDs in the clinical data
    clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
    clinicalInd <- match(completeIDs, clinicalIDs)
    clinicalData <- clinicalData[clinicalInd, ]
    # Subset clinical data and completeData by clinicalData
    clinicalData <- clinicalData[rowSums(is.na(clinicalData)) < ncol(clinicalData),]
    clinicalIDs <- as.character(clinicalData$bcr_patient_barcode)
    completeInd <- match(clinicalIDs, completeIDs)
    
    ###### Run combat function 
    for (i in 1:numViews) {
      temp.completeData <- completeData[[i]]
      temp.completeData <- temp.completeData[,completeInd]
      temp.2.completeData <- temp.completeData[rowSums(temp.completeData) != 0,]
      temp.modcombat <- model.matrix(~1, data = clinicalData)
      temp.batch <- clinicalData$gender
      temp_combat = ComBat(dat=temp.2.completeData, batch=temp.batch, mod=temp.modcombat, par.prior=TRUE, prior.plots=FALSE)
      completeData[[i]] <- temp_combat
    }
    
    
    completeIds <- transformUnionId(colnames(completeData[[1]]) )
    
    return(list(first = completeIds, second = unionIds))
    
  } else { 
    
    dataTypes <- c("methyl", "mirna", "mrna")
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
    
    
    
    
    # extract all cases which appear in all of the data types (intersection)
    completeData <- columnIntersection(cases) # now ncol in each matrix is the same and identical to each other. 
    
    # subset the clinical data so that it corresponds to individuals in the complete data
    completeIds <- colnames(completeData[[1]])   
    
    # Extract all cases which appear in at least one of the data types
    unionData <- columnUnion(cases)
    
    # Subset the clinical data so that it corresponds to individuals
    # in the union data
    unionIds <- colnames(unionData[[1]])
    
    return(list(first = completeIds, second = unionIds))
    
  }
  
}

# Lod cancer ids
brca_ids <-loadIDs(cancer = 'BRCA')
kirc_ids <-loadIDs(cancer = 'KIRC')
lihc_ids <-loadIDs(cancer = 'LIHC')
luad_ids <-loadIDs(cancer = 'LUAD')

##########################################################################################################
# BRCA, compare clusters from complete, compete imputed, and union

# get brca IDs
brca_complete_ids <- brca_ids[[1]]
brca_union_ids <- brca_ids[[2]]

# Get labels for each clustering type, from the intersection for brca
brca_com5_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_1_1.txt'))))
brca_com5_iclust <- as.factor(t(read.table(paste0(completeFolder, '/1_1_2.txt'))))
brca_com5_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_1_3.txt'))))

brca_com4_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_2_1.txt'))))
brca_com4_iclust <- as.factor(t(read.table(paste0(completeFolder, '/1_2_2.txt'))))
brca_com4_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_2_3.txt'))))

brca_com3_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_3_1.txt'))))
brca_com3_iclust <- as.factor(t(read.table(paste0(completeFolder, '/1_3_2.txt'))))
brca_com3_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_3_3.txt'))))

brca_com2_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_4_1.txt'))))
brca_com2_iclust <- as.factor(t(read.table(paste0(completeFolder, '/1_4_2.txt'))))
brca_com2_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_4_3.txt'))))
############################################################################################3
# get labels for brca intersection imputed
# knn cluster 1 - 4
brca_cluster5_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster1_knn_hier.txt'))))
brca_cluster4_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster2_knn_hier.txt'))))
brca_cluster3_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster3_knn_hier.txt'))))
brca_cluster2_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster4_knn_hier.txt'))))

brca_cluster5_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster1_knn_iclust.txt'))))
brca_cluster4_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster2_knn_iclust.txt'))))
brca_cluster3_com_knn_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster3_knn_iclust.txt'))))
brca_cluster2_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster4_knn_iclust.txt'))))

brca_cluster5_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster1_knn_snf.txt'))))
brca_cluster4_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster2_knn_snf.txt'))))
brca_cluster3_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster3_knn_snf.txt'))))
brca_cluster2_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
brca_cluster5_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster1_lls_hier.txt'))))
brca_cluster4_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster2_lls_hier.txt'))))
brca_cluster3_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster3_lls_hier.txt'))))
brca_cluster2_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster4_lls_hier.txt'))))

brca_cluster5_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster1_lls_iclust.txt'))))
brca_cluster4_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster2_lls_iclust.txt'))))
brca_cluster4_com_lls_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster3_lls_iclust.txt'))))
brca_cluster2_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster4_lls_iclust.txt'))))

brca_cluster5_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster1_lls_snf.txt'))))
brca_cluster4_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster2_lls_snf.txt'))))
brca_cluster3_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster3_lls_snf.txt'))))
brca_cluster2_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
brca_cluster5_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster1_lsa_hier.txt'))))
brca_cluster4_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster2_lsa_hier.txt'))))
brca_cluster3_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster3_lsa_hier.txt'))))
brca_cluster2_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster4_lsa_hier.txt'))))

brca_cluster5_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster1_lsa_iclust.txt'))))
brca_cluster4_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster2_lsa_iclust.txt'))))
brca_cluster3_com_lsa_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster3_lsa_iclust.txt'))))
brca_cluster2_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster4_lsa_iclust.txt'))))

brca_cluster5_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster1_lsa_snf.txt'))))
brca_cluster2_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster2_lsa_snf.txt'))))
brca_cluster3_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster3_lsa_snf.txt'))))
brca_cluster2_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
brca_cluster5_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster1_rand_hier.txt'))))
brca_cluster4_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster2_rand_hier.txt'))))
brca_cluster3_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster3_rand_hier.txt'))))
brca_cluster2_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster4_rand_hier.txt'))))

brca_cluster5_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster1_rand_iclust.txt'))))
brca_cluster4_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster2_rand_iclust.txt'))))
brca_cluster3_com_rand_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster3_rand_iclust.txt'))))
brca_cluster2_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster4_rand_iclust.txt'))))

brca_cluster5_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster1_rand_snf.txt'))))
brca_cluster4_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster2_rand_snf.txt'))))
brca_cluster3_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster3_rand_snf.txt'))))
brca_cluster2_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/rand', '/cluster4_rand_snf.txt'))))


# similarity
brca_cluster5_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster1_med_snf.txt'))))
brca_cluster4_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster2_med_snf.txt'))))
brca_cluster3_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster3_med_snf.txt'))))
brca_cluster2_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster4_med_snf.txt'))))


brca_cluster5_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster1_regress_snf.txt'))))
brca_cluster4_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster2_regress_snf.txt'))))
brca_cluster3_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster3_regress_snf.txt'))))
brca_cluster2_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster4_regress_snf.txt'))))

brca_cluster5_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster1_self_snf.txt'))))
brca_cluster4_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster2_self_snf.txt'))))
brca_cluster3_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster3_self_snf.txt'))))
brca_cluster2_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/cluster4_self_snf.txt'))))

#########################################################################################################
# Get labels for union

# knn cluster 1 - 4
brca_cluster5_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster1_knn_hier.txt'))))
brca_cluster4_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster2_knn_hier.txt'))))
brca_cluster3_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster3_knn_hier.txt'))))
brca_cluster2_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster4_knn_hier.txt'))))

brca_cluster5_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster1_knn_iclust.txt'))))
brca_cluster4_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster2_knn_iclust.txt'))))
brca_cluster3_union_knn_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster3_knn_iclust.txt'))))
brca_cluster2_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster4_knn_iclust.txt'))))

brca_cluster5_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster1_knn_snf.txt'))))
brca_cluster4_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster2_knn_snf.txt'))))
brca_cluster3_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster3_knn_snf.txt'))))
brca_cluster2_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
brca_cluster5_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster1_lls_hier.txt'))))
brca_cluster4_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster2_lls_hier.txt'))))
brca_cluster3_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster3_lls_hier.txt'))))
brca_cluster2_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster4_lls_hier.txt'))))

brca_cluster5_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster1_lls_iclust.txt'))))
brca_cluster4_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster2_lls_iclust.txt'))))
brca_cluster3_union_lls_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster3_lls_iclust.txt'))))
brca_cluster2_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster4_lls_iclust.txt'))))

brca_cluster5_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster1_lls_snf.txt'))))
brca_cluster4_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster2_lls_snf.txt'))))
brca_cluster3_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster3_lls_snf.txt'))))
brca_cluster2_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
brca_cluster5_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster1_lsa_hier.txt'))))
brca_cluster4_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster2_lsa_hier.txt'))))
brca_cluster3_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster3_lsa_hier.txt'))))
brca_cluster2_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster4_lsa_hier.txt'))))

brca_cluster5_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster1_lsa_iclust.txt'))))
brca_cluster4_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster2_lsa_iclust.txt'))))
brca_cluster3_union_lsa_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster3_lsa_iclust.txt'))))
brca_cluster2_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster4_lsa_iclust.txt'))))

brca_cluster5_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster1_lsa_snf.txt'))))
brca_cluster4_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster2_lsa_snf.txt'))))
brca_cluster3_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster3_lsa_snf.txt'))))
brca_cluster2_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
brca_cluster5_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster1_rand_hier.txt'))))
brca_cluster4_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster2_rand_hier.txt'))))
brca_cluster3_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster3_rand_hier.txt'))))
brca_cluster2_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster4_rand_hier.txt'))))

brca_cluster5_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster1_rand_iclust.txt'))))
brca_cluster4_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster2_rand_iclust.txt'))))
brca_cluster3_union_rand_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster3_rand_iclust.txt'))))
brca_cluster2_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster4_rand_iclust.txt'))))

brca_cluster5_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster1_rand_snf.txt'))))
brca_cluster4_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster2_rand_snf.txt'))))
brca_cluster3_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster3_rand_snf.txt'))))
brca_cluster2_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/rand', '/cluster4_rand_snf.txt'))))


# similarity
brca_cluster5_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster1_med_snf.txt'))))
brca_cluster4_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster2_med_snf.txt'))))
brca_cluster3_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster3_med_snf.txt'))))
brca_cluster2_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster4_med_snf.txt'))))


brca_cluster5_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster1_regress_snf.txt'))))
brca_cluster4_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster2_regress_snf.txt'))))
brca_cluster3_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster3_regress_snf.txt'))))
brca_cluster2_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster4_regress_snf.txt'))))

brca_cluster5_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster1_self_snf.txt'))))
brca_cluster4_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster2_self_snf.txt'))))
brca_cluster3_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster3_self_snf.txt'))))
brca_cluster2_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/cluster4_self_snf.txt'))))


##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical 
brca_hier_stat <- rbind(
  append('brca_com_hier_5', summary(as.factor(brca_com5_hier))),
  append('brca_com_hier_4', summary(as.factor(brca_com4_hier))),
  append('brca_com_hier_3', summary(as.factor(brca_com3_hier))),
  append('brca_com_hier_2', summary(as.factor(brca_com2_hier)),
  append('brca_com_hier_rand_5', summary(as.factor(brca_cluster5_com_rand_hier))),
  append('brca_union_hier_knn', summary(as.factor(brca_union_hier_knn))),
  append('brca_union_hier_lls', summary(as.factor(brca_union_hier_lls))),
  append('brca_union_hier_lsa', summary(as.factor(brca_union_hier_lsa))),
  append('brca_union_hier_rand', summary(as.factor(brca_union_hier_rand)))
  
)

# Icluster
brca_iclust_stat <- rbind(
  append('brca_com_iclust', summary(as.factor(brca_com_iclus))),
  append('brca_com_iclust_knn', summary(as.factor(brca_com_iclus_knn))),
  append('brca_com_iclust_lls', summary(as.factor(brca_com_iclus_lls))),
  append('brca_com_iclust_lsa', summary(as.factor(brca_com_iclus_lsa))),
  append('brca_com_iclust_rand', summary(as.factor(brca_com_iclus_rand))),
  append('brca_union_iclust_knn', summary(as.factor(brca_union_iclus_knn))),
  append('brca_union_iclust_lls', summary(as.factor(brca_union_iclus_lls))),
  append('brca_union_iclust_lsa', summary(as.factor(brca_union_iclus_lsa))),
  append('brca_union_iclust_rand', summary(as.factor(brca_union_iclus_rand)))
  
)


# SNF
brca_snf_stat <- rbind(
  append('brca_com_snf', summary(as.factor(brca_com_snf))),
  append('brca_com_snf_self', summary(as.factor(brca_com_snf_self))),
  append('brca_com_snf_med', summary(as.factor(brca_com_snf_med))),
  append('brca_com_snf_reg', summary(as.factor(brca_com_snf_reg))),
  append('brca_union_snf_self', summary(as.factor(brca_union_snf_self))),
  append('brca_union_snf_med', summary(as.factor(brca_union_snf_med))),
  append('brca_union_snf_reg', summary(as.factor(brca_union_snf_reg)))
)



# Merge ids with labels and get clinical information for each cluster

# look at groups of ids in complete, complete_impute, and union
# join ids and then join each label groupd
brca_com_hier <- as.data.frame(cbind(label = brca_com_hier, id = brca_complete_ids))
brca_com_hier_knn <- as.data.frame(cbind(label = brca_com_hier_knn, id = brca_complete_ids))
brca_com_hier_lls <- as.data.frame(cbind(label = brca_com_hier_lls, id = brca_complete_ids))
brca_com_hier_lsa <- as.data.frame(cbind(label = brca_com_hier_lsa, id = brca_complete_ids))
brca_com_hier_rand <-  as.data.frame(cbind(label = brca_union_hier_rand, id = brca_union_ids))
brca_union_hier_knn <- as.data.frame(cbind(label = brca_union_hier_knn, id = brca_union_ids))
brca_union_hier_lls <- as.data.frame(cbind(label = brca_union_hier_lls, id = brca_union_ids))
brca_union_hier_lsa <- as.data.frame(cbind(label = brca_union_hier_lsa, id = brca_union_ids))
brca_union_hier_rand <-  as.data.frame(cbind(label = brca_union_hier_rand, id = brca_union_ids))

brca_com_iclus <- as.data.frame(cbind(label = brca_com_iclus, id = brca_complete_ids))
brca_com_iclus_knn <- as.data.frame(cbind(label = brca_com_iclus_knn, id = brca_complete_ids))
brca_com_iclus_lls <- as.data.frame(cbind(label = brca_com_iclus_lls, id = brca_complete_ids))
brca_com_iclus_lsa <- as.data.frame(cbind(label = brca_com_iclus_lsa, id = brca_complete_ids))
brca_com_iclus_rand <-  as.data.frame(cbind(label = brca_union_iclus_rand, id = brca_union_ids))
brca_union_iclus_knn <- as.data.frame(cbind(label = brca_union_iclus_knn, id = brca_union_ids))
brca_union_iclus_lls <- as.data.frame(cbind(label = brca_union_iclus_lls, id = brca_union_ids))
brca_union_iclus_lsa <- as.data.frame(cbind(label = brca_union_iclus_lsa, id = brca_union_ids))
brca_union_iclus_rand <-  as.data.frame(cbind(label = brca_union_iclus_rand, id = brca_union_ids))

brca_com_snf <- as.data.frame(cbind(label = brca_com_snf, id = brca_complete_ids))
brca_com_snf_knn <- as.data.frame(cbind(label = brca_com_snf_knn, id = brca_complete_ids))
brca_com_snf_lls <- as.data.frame(cbind(label = brca_com_snf_lls, id = brca_complete_ids))
brca_com_snf_lsa <- as.data.frame(cbind(label = brca_com_snf_lsa, id = brca_complete_ids))
brca_com_snf_rand <-  as.data.frame(cbind(label = brca_union_snf_rand, id = brca_union_ids))
brca_com_snf_self <- as.data.frame(cbind(label = brca_com_snf_self, id = brca_complete_ids))
brca_com_snf_med <- as.data.frame(cbind(label = brca_com_snf_med, id = brca_complete_ids))
brca_com_snf_reg <- as.data.frame(cbind(label = brca_com_snf_reg, id = brca_complete_ids))
brca_union_snf_knn <- as.data.frame(cbind(label = brca_union_snf_knn, id = brca_union_ids))
brca_union_snf_lls <- as.data.frame(cbind(label = brca_union_snf_lls, id = brca_union_ids))
brca_union_snf_lsa <- as.data.frame(cbind(label = brca_union_snf_lsa, id = brca_union_ids))
brca_union_snf_rand <-  as.data.frame(cbind(label = brca_union_snf_rand, id = brca_union_ids))
brca_union_snf_self <-  as.data.frame(cbind(label = brca_union_snf_self, id = brca_union_ids))
brca_union_snf_med <-  as.data.frame(cbind(label = brca_union_snf_med, id = brca_union_ids))
brca_union_snf_reg <-  as.data.frame(cbind(label = brca_union_snf_reg, id = brca_union_ids))

findCommonCluster <- function(data1, data2) {
  numClus <- 5
  data_intersect_union <- data2[data2$id %in% data1$id,] # intersection of union
  data_intersect <- left_join(data_intersect_union, data1, by = 'id') # join intersection of union with intersection
  common_clusters <- matrix(,0,5)
  cluster_compare <- matrix(,0,5)
  for(i in 1:numClus){
    sub_intersect <- data_intersect[data_intersect$label.y == i,]
    for(j in 1:numClus){
      sub_union <- data_intersect[data_intersect$label.x == j,]
      cluster_compare[j] <- round((sum(sub_intersect$id %in% sub_union$id)/nrow(sub_intersect))*100, 2)
    }
    common_clusters <- rbind(common_clusters, cluster_compare)
  }
  return(common_clusters)
}

brca_hier_knn_common <- findCommonCluster(brca_com_hier_knn, brca_union_hier_knn)
brca_hier_lls_common <- findCommonCluster(brca_com_hier_lls, brca_union_hier_lls)
brca_hier_lsa_common <- findCommonCluster(brca_com_hier_lsa, brca_union_hier_lsa)
brca_hier_rand_common <- findCommonCluster(brca_com_hier_rand, brca_union_hier_rand)

brca_iclus_knn_common <- findCommonCluster(brca_com_iclus_knn, brca_union_iclus_knn)
brca_iclus_lls_common <- findCommonCluster(brca_com_iclus_lls, brca_union_iclus_lls)
brca_iclus_lsa_common <- findCommonCluster(brca_com_iclus_lsa, brca_union_iclus_lsa)
brca_iclus_rand_common <- findCommonCluster(brca_com_iclus_rand, brca_union_iclus_rand)

brca_snf_knn_common <- findCommonCluster(brca_com_snf_knn, brca_union_snf_knn)
brca_snf_lls_common <- findCommonCluster(brca_com_snf_lls, brca_union_snf_lls)
brca_snf_lsa_common <- findCommonCluster(brca_com_snf_lsa, brca_union_snf_lsa)
brca_snf_rand_common <- findCommonCluster(brca_com_snf_rand, brca_union_snf_rand)
brca_snf_self_common <- findCommonCluster(brca_com_snf_self, brca_union_snf_self)
brca_snf_reg_common <- findCommonCluster(brca_com_snf_reg, brca_union_snf_reg)
brca_snf_med_common <- findCommonCluster(brca_com_snf_med, brca_union_snf_med)
