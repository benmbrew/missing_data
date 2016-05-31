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
lusc_ids <-loadIDs(cancer = 'LUSC')


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
brca_cluster3_com_lls_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/brca/lls', '/cluster3_lls_iclust.txt'))))
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
brca_cluster4_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/brca/lsa', '/cluster2_lsa_snf.txt'))))
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

# Hierarchical############################################################################################

# Hierarchical 5 clusters
brca_hier_stat_5 <- rbind(
  append('brca_com_hier_5', summary(as.factor(brca_com5_hier))),
  append('brca_com_hier_knn_5', summary(as.factor(brca_cluster5_com_knn_hier))),
  append('brca_com_hier_lls_5', summary(as.factor(brca_cluster5_com_lls_hier))),
  append('brca_com_hier_lsa_5', summary(as.factor(brca_cluster5_com_lsa_hier))),
  append('brca_com_hier_rand_5', summary(as.factor(brca_cluster5_com_rand_hier))),
  append('brca_union_hier_knn_5', summary(as.factor(brca_cluster5_union_knn_hier))),
  append('brca_union_hier_lls_5', summary(as.factor(brca_cluster5_union_lls_hier))),
  append('brca_union_hier_lsa_5', summary(as.factor(brca_cluster5_union_lsa_hier))),
  append('brca_union_hier_rand_5', summary(as.factor(brca_cluster5_union_rand_hier)))
)

# Hierarchical 4 clusters
brca_hier_stat_4 <- rbind(
  append('brca_com_hier_4', summary(as.factor(brca_com4_hier))),
  append('brca_com_hier_knn_4', summary(as.factor(brca_cluster4_com_knn_hier))),
  append('brca_com_hier_lls_4', summary(as.factor(brca_cluster4_com_lls_hier))),
  append('brca_com_hier_lsa_4', summary(as.factor(brca_cluster4_com_lsa_hier))),
  append('brca_com_hier_rand_4', summary(as.factor(brca_cluster4_com_rand_hier))),
  append('brca_union_hier_knn_4', summary(as.factor(brca_cluster4_union_knn_hier))),
  append('brca_union_hier_lls_4', summary(as.factor(brca_cluster4_union_lls_hier))),
  append('brca_union_hier_lsa_4', summary(as.factor(brca_cluster4_union_lsa_hier))),
  append('brca_union_hier_rand_4', summary(as.factor(brca_cluster4_union_rand_hier)))
)

# Hierarchical 3 clusters
brca_hier_stat_3 <- rbind(
  append('brca_com_hier_3', summary(as.factor(brca_com3_hier))),
  append('brca_com_hier_knn_3', summary(as.factor(brca_cluster3_com_knn_hier))),
  append('brca_com_hier_lls_3', summary(as.factor(brca_cluster3_com_lls_hier))),
  append('brca_com_hier_lsa_3', summary(as.factor(brca_cluster3_com_lsa_hier))),
  append('brca_com_hier_rand_3', summary(as.factor(brca_cluster3_com_rand_hier))),
  append('brca_union_hier_knn_3', summary(as.factor(brca_cluster3_union_knn_hier))),
  append('brca_union_hier_lls_3', summary(as.factor(brca_cluster3_union_lls_hier))),
  append('brca_union_hier_lsa_3', summary(as.factor(brca_cluster3_union_lsa_hier))),
  append('brca_union_hier_rand_3', summary(as.factor(brca_cluster3_union_rand_hier)))
)

# Hierarchical 2 clusters
brca_hier_stat_2 <- rbind(
  append('brca_com_hier_2', summary(as.factor(brca_com2_hier))),
  append('brca_com_hier_knn_2', summary(as.factor(brca_cluster2_com_knn_hier))),
  append('brca_com_hier_lls_2', summary(as.factor(brca_cluster2_com_lls_hier))),
  append('brca_com_hier_lsa_2', summary(as.factor(brca_cluster2_com_lsa_hier))),
  append('brca_com_hier_rand_2', summary(as.factor(brca_cluster2_com_rand_hier))),
  append('brca_union_hier_knn_2', summary(as.factor(brca_cluster2_union_knn_hier))),
  append('brca_union_hier_lls_2', summary(as.factor(brca_cluster2_union_lls_hier))),
  append('brca_union_hier_lsa_2', summary(as.factor(brca_cluster2_union_lsa_hier))),
  append('brca_union_hier_rand_2', summary(as.factor(brca_cluster2_union_rand_hier)))
)


# icluster############################################################################################

# icluster 5 clusters
brca_iclust_stat_5 <- rbind(
  append('brca_com_iclust_5', summary(as.factor(brca_com5_iclust))),
  append('brca_com_iclust_knn_5', summary(as.factor(brca_cluster5_com_knn_iclust))),
  append('brca_com_iclust_lls_5', summary(as.factor(brca_cluster5_com_lls_iclust))),
  append('brca_com_iclust_lsa_5', summary(as.factor(brca_cluster5_com_lsa_iclust))),
  append('brca_com_iclust_rand_5', summary(as.factor(brca_cluster5_com_rand_iclust))),
  append('brca_union_iclust_knn_5', summary(as.factor(brca_cluster5_union_knn_iclust))),
  append('brca_union_iclust_lls_5', summary(as.factor(brca_cluster5_union_lls_iclust))),
  append('brca_union_iclust_lsa_5', summary(as.factor(brca_cluster5_union_lsa_iclust))),
  append('brca_union_iclust_rand_5', summary(as.factor(brca_cluster5_union_rand_iclust)))
)

# icluster 4 clusters
brca_iclust_stat_4 <- rbind(
  append('brca_com_iclust_4', summary(as.factor(brca_com4_iclust))),
  append('brca_com_iclust_knn_4', summary(as.factor(brca_cluster4_com_knn_iclust))),
  append('brca_com_iclust_lls_4', summary(as.factor(brca_cluster4_com_lls_iclust))),
  append('brca_com_iclust_lsa_4', summary(as.factor(brca_cluster4_com_lsa_iclust))),
  append('brca_com_iclust_rand_4', summary(as.factor(brca_cluster4_com_rand_iclust))),
  append('brca_union_iclust_knn_4', summary(as.factor(brca_cluster4_union_knn_iclust))),
  append('brca_union_iclust_lls_4', summary(as.factor(brca_cluster4_union_lls_iclust))),
  append('brca_union_iclust_lsa_4', summary(as.factor(brca_cluster4_union_lsa_iclust))),
  append('brca_union_iclust_rand_4', summary(as.factor(brca_cluster4_union_rand_iclust)))
)

# icluster 3 clusters
brca_iclust_stat_3 <- rbind(
  append('brca_com_iclust_3', summary(as.factor(brca_com3_iclust))),
  append('brca_com_iclust_knn_3', summary(as.factor(brca_cluster3_com_knn_iclust))),
  append('brca_com_iclust_lls_3', summary(as.factor(brca_cluster3_com_lls_iclust))),
  append('brca_com_iclust_lsa_3', summary(as.factor(brca_cluster3_com_lsa_iclust))),
  append('brca_com_iclust_rand_3', summary(as.factor(brca_cluster3_com_rand_iclust))),
  append('brca_union_iclust_knn_3', summary(as.factor(brca_cluster3_union_knn_iclust))),
  append('brca_union_iclust_lls_3', summary(as.factor(brca_cluster3_union_lls_iclust))),
  append('brca_union_iclust_lsa_3', summary(as.factor(brca_cluster3_union_lsa_iclust))),
  append('brca_union_iclust_rand_3', summary(as.factor(brca_cluster3_union_rand_iclust)))
)

# icluster 2 clusters
brca_iclust_stat_2 <- rbind(
  append('brca_com_iclust_2', summary(as.factor(brca_com2_iclust))),
  append('brca_com_iclust_knn_2', summary(as.factor(brca_cluster2_com_knn_iclust))),
  append('brca_com_iclust_lls_2', summary(as.factor(brca_cluster2_com_lls_iclust))),
  append('brca_com_iclust_lsa_2', summary(as.factor(brca_cluster2_com_lsa_iclust))),
  append('brca_com_iclust_rand_2', summary(as.factor(brca_cluster2_com_rand_iclust))),
  append('brca_union_iclust_knn_2', summary(as.factor(brca_cluster2_union_knn_iclust))),
  append('brca_union_iclust_lls_2', summary(as.factor(brca_cluster2_union_lls_iclust))),
  append('brca_union_iclust_lsa_2', summary(as.factor(brca_cluster2_union_lsa_iclust))),
  append('brca_union_iclust_rand_2', summary(as.factor(brca_cluster2_union_rand_iclust)))
)


# snf############################################################################################

# snf 5 clusters
brca_snf_stat_5 <- rbind(
  
  append('brca_com_snf_5', summary(as.factor(brca_com5_snf))),
  append('brca_com_snf_knn_5', summary(as.factor(brca_cluster5_com_knn_snf))),
  append('brca_com_snf_lls_5', summary(as.factor(brca_cluster5_com_lls_snf))),
  append('brca_com_snf_lsa_5', summary(as.factor(brca_cluster5_com_lsa_snf))),
  append('brca_com_snf_rand_5', summary(as.factor(brca_cluster5_com_rand_snf))),
  append('brca_com_snf_reg_5', summary(as.factor(brca_cluster5_com_reg_snf))),
  append('brca_com_snf_med_5', summary(as.factor(brca_cluster5_com_med_snf))),
  append('brca_com_snf_self_5', summary(as.factor(brca_cluster5_com_self_snf))),
  append('brca_union_snf_knn_5', summary(as.factor(brca_cluster5_union_knn_snf))),
  append('brca_union_snf_lls_5', summary(as.factor(brca_cluster5_union_lls_snf))),
  append('brca_union_snf_lsa_5', summary(as.factor(brca_cluster5_union_lsa_snf))),
  append('brca_union_snf_rand_5', summary(as.factor(brca_cluster5_union_rand_snf))),
  append('brca_union_snf_reg_5', summary(as.factor(brca_cluster5_union_reg_snf))),
  append('brca_union_snf_med_5', summary(as.factor(brca_cluster5_union_med_snf))),
  append('brca_union_snf_self_5', summary(as.factor(brca_cluster5_union_self_snf)))
  
)

# snf 4 clusters
brca_snf_stat_4 <- rbind(
  
  append('brca_com_snf_4', summary(as.factor(brca_com4_snf))),
  append('brca_com_snf_knn_4', summary(as.factor(brca_cluster4_com_knn_snf))),
  append('brca_com_snf_lls_4', summary(as.factor(brca_cluster4_com_lls_snf))),
  append('brca_com_snf_lsa_4', summary(as.factor(brca_cluster4_com_lsa_snf))),
  append('brca_com_snf_rand_4', summary(as.factor(brca_cluster4_com_rand_snf))),
  append('brca_com_snf_reg_4', summary(as.factor(brca_cluster4_com_reg_snf))),
  append('brca_com_snf_med_4', summary(as.factor(brca_cluster4_com_med_snf))),
  append('brca_com_snf_self_4', summary(as.factor(brca_cluster4_com_self_snf))),
  append('brca_union_snf_knn_4', summary(as.factor(brca_cluster4_union_knn_snf))),
  append('brca_union_snf_lls_4', summary(as.factor(brca_cluster4_union_lls_snf))),
  append('brca_union_snf_lsa_4', summary(as.factor(brca_cluster4_union_lsa_snf))),
  append('brca_union_snf_rand_4', summary(as.factor(brca_cluster4_union_rand_snf))),
  append('brca_union_snf_reg_4', summary(as.factor(brca_cluster4_union_reg_snf))),
  append('brca_union_snf_med_4', summary(as.factor(brca_cluster4_union_med_snf))),
  append('brca_union_snf_self_4', summary(as.factor(brca_cluster4_union_self_snf)))
  
)


# snf 3 clusters
brca_snf_stat_3 <- rbind(
  
  append('brca_com_snf_3', summary(as.factor(brca_com3_snf))),
  append('brca_com_snf_knn_3', summary(as.factor(brca_cluster3_com_knn_snf))),
  append('brca_com_snf_lls_3', summary(as.factor(brca_cluster3_com_lls_snf))),
  append('brca_com_snf_lsa_3', summary(as.factor(brca_cluster3_com_lsa_snf))),
  append('brca_com_snf_rand_3', summary(as.factor(brca_cluster3_com_rand_snf))),
  append('brca_com_snf_reg_3', summary(as.factor(brca_cluster3_com_reg_snf))),
  append('brca_com_snf_med_3', summary(as.factor(brca_cluster3_com_med_snf))),
  append('brca_com_snf_self_3', summary(as.factor(brca_cluster3_com_self_snf))),
  append('brca_union_snf_knn_3', summary(as.factor(brca_cluster3_union_knn_snf))),
  append('brca_union_snf_lls_3', summary(as.factor(brca_cluster3_union_lls_snf))),
  append('brca_union_snf_lsa_3', summary(as.factor(brca_cluster3_union_lsa_snf))),
  append('brca_union_snf_rand_3', summary(as.factor(brca_cluster3_union_rand_snf))),
  append('brca_union_snf_reg_3', summary(as.factor(brca_cluster3_union_reg_snf))),
  append('brca_union_snf_med_3', summary(as.factor(brca_cluster3_union_med_snf))),
  append('brca_union_snf_self_3', summary(as.factor(brca_cluster3_union_self_snf)))
  
)

# snf 2 clusters
brca_snf_stat_2 <- rbind(
  
  append('brca_com_snf_2', summary(as.factor(brca_com2_snf))),
  append('brca_com_snf_knn_2', summary(as.factor(brca_cluster2_com_knn_snf))),
  append('brca_com_snf_lls_2', summary(as.factor(brca_cluster2_com_lls_snf))),
  append('brca_com_snf_lsa_2', summary(as.factor(brca_cluster2_com_lsa_snf))),
  append('brca_com_snf_rand_2', summary(as.factor(brca_cluster2_com_rand_snf))),
  append('brca_com_snf_reg_2', summary(as.factor(brca_cluster2_com_reg_snf))),
  append('brca_com_snf_med_2', summary(as.factor(brca_cluster2_com_med_snf))),
  append('brca_com_snf_self_2', summary(as.factor(brca_cluster2_com_self_snf))),
  append('brca_union_snf_knn_2', summary(as.factor(brca_cluster2_union_knn_snf))),
  append('brca_union_snf_lls_2', summary(as.factor(brca_cluster2_union_lls_snf))),
  append('brca_union_snf_lsa_2', summary(as.factor(brca_cluster2_union_lsa_snf))),
  append('brca_union_snf_rand_2', summary(as.factor(brca_cluster2_union_rand_snf))),
  append('brca_union_snf_reg_2', summary(as.factor(brca_cluster2_union_reg_snf))),
  append('brca_union_snf_med_2', summary(as.factor(brca_cluster2_union_med_snf))),
  append('brca_union_snf_self_2', summary(as.factor(brca_cluster2_union_self_snf)))
  
)

#################################################################################################

##########################################################################################################
# kirc, compare clusters from complete, compete imputed, and union

# get kirc IDs
kirc_complete_ids <- kirc_ids[[1]]
kirc_union_ids <- kirc_ids[[2]]

# Get labels for each clustering type, from the intersection for kirc
kirc_com5_hier <- as.factor(t(read.table(paste0(completeFolder, '/2_1_1.txt'))))
kirc_com5_iclust <- as.factor(t(read.table(paste0(completeFolder, '/2_1_2.txt'))))
kirc_com5_snf <- as.factor(t(read.table(paste0(completeFolder, '/2_1_3.txt'))))

kirc_com4_hier <- as.factor(t(read.table(paste0(completeFolder, '/2_2_1.txt'))))
kirc_com4_iclust <- as.factor(t(read.table(paste0(completeFolder, '/2_2_2.txt'))))
kirc_com4_snf <- as.factor(t(read.table(paste0(completeFolder, '/2_2_3.txt'))))

kirc_com3_hier <- as.factor(t(read.table(paste0(completeFolder, '/2_3_1.txt'))))
kirc_com3_iclust <- as.factor(t(read.table(paste0(completeFolder, '/2_3_2.txt'))))
kirc_com3_snf <- as.factor(t(read.table(paste0(completeFolder, '/2_3_3.txt'))))

kirc_com2_hier <- as.factor(t(read.table(paste0(completeFolder, '/2_4_1.txt'))))
kirc_com2_iclust <- as.factor(t(read.table(paste0(completeFolder, '/2_4_2.txt'))))
kirc_com2_snf <- as.factor(t(read.table(paste0(completeFolder, '/2_4_3.txt'))))
############################################################################################3
# get labels for kirc intersection imputed
# knn cluster 1 - 4
kirc_cluster5_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster1_knn_hier.txt'))))
kirc_cluster4_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster2_knn_hier.txt'))))
kirc_cluster3_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster3_knn_hier.txt'))))
kirc_cluster2_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster4_knn_hier.txt'))))

kirc_cluster5_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster1_knn_iclust.txt'))))
kirc_cluster4_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster2_knn_iclust.txt'))))
kirc_cluster3_com_knn_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster3_knn_iclust.txt'))))
kirc_cluster2_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster4_knn_iclust.txt'))))

kirc_cluster5_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster1_knn_snf.txt'))))
kirc_cluster4_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster2_knn_snf.txt'))))
kirc_cluster3_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster3_knn_snf.txt'))))
kirc_cluster2_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
kirc_cluster5_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster1_lls_hier.txt'))))
kirc_cluster4_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster2_lls_hier.txt'))))
kirc_cluster3_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster3_lls_hier.txt'))))
kirc_cluster2_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster4_lls_hier.txt'))))

kirc_cluster5_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster1_lls_iclust.txt'))))
kirc_cluster4_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster2_lls_iclust.txt'))))
kirc_cluster3_com_lls_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster3_lls_iclust.txt'))))
kirc_cluster2_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster4_lls_iclust.txt'))))

kirc_cluster5_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster1_lls_snf.txt'))))
kirc_cluster4_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster2_lls_snf.txt'))))
kirc_cluster3_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster3_lls_snf.txt'))))
kirc_cluster2_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
kirc_cluster5_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster1_lsa_hier.txt'))))
kirc_cluster4_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster2_lsa_hier.txt'))))
kirc_cluster3_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster3_lsa_hier.txt'))))
kirc_cluster2_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster4_lsa_hier.txt'))))

kirc_cluster5_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster1_lsa_iclust.txt'))))
kirc_cluster4_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster2_lsa_iclust.txt'))))
kirc_cluster3_com_lsa_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster3_lsa_iclust.txt'))))
kirc_cluster2_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster4_lsa_iclust.txt'))))

kirc_cluster5_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster1_lsa_snf.txt'))))
kirc_cluster4_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster2_lsa_snf.txt'))))
kirc_cluster3_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster3_lsa_snf.txt'))))
kirc_cluster2_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
kirc_cluster5_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster1_rand_hier.txt'))))
kirc_cluster4_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster2_rand_hier.txt'))))
kirc_cluster3_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster3_rand_hier.txt'))))
kirc_cluster2_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster4_rand_hier.txt'))))

kirc_cluster5_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster1_rand_iclust.txt'))))
kirc_cluster4_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster2_rand_iclust.txt'))))
kirc_cluster3_com_rand_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster3_rand_iclust.txt'))))
kirc_cluster2_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster4_rand_iclust.txt'))))

kirc_cluster5_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster1_rand_snf.txt'))))
kirc_cluster4_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster2_rand_snf.txt'))))
kirc_cluster3_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster3_rand_snf.txt'))))
kirc_cluster2_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/kirc/rand', '/cluster4_rand_snf.txt'))))


# similarity
kirc_cluster5_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster1_med_snf.txt'))))
kirc_cluster4_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster2_med_snf.txt'))))
kirc_cluster3_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster3_med_snf.txt'))))
kirc_cluster2_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster4_med_snf.txt'))))


kirc_cluster5_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster1_regress_snf.txt'))))
kirc_cluster4_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster2_regress_snf.txt'))))
kirc_cluster3_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster3_regress_snf.txt'))))
kirc_cluster2_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster4_regress_snf.txt'))))

kirc_cluster5_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster1_self_snf.txt'))))
kirc_cluster4_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster2_self_snf.txt'))))
kirc_cluster3_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster3_self_snf.txt'))))
kirc_cluster2_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/cluster4_self_snf.txt'))))

#########################################################################################################
# Get labels for union

# knn cluster 1 - 4
kirc_cluster5_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster1_knn_hier.txt'))))
kirc_cluster4_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster2_knn_hier.txt'))))
kirc_cluster3_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster3_knn_hier.txt'))))
kirc_cluster2_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster4_knn_hier.txt'))))

kirc_cluster5_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster1_knn_iclust.txt'))))
kirc_cluster4_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster2_knn_iclust.txt'))))
kirc_cluster3_union_knn_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster3_knn_iclust.txt'))))
kirc_cluster2_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster4_knn_iclust.txt'))))

kirc_cluster5_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster1_knn_snf.txt'))))
kirc_cluster4_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster2_knn_snf.txt'))))
kirc_cluster3_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster3_knn_snf.txt'))))
kirc_cluster2_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
kirc_cluster5_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster1_lls_hier.txt'))))
kirc_cluster4_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster2_lls_hier.txt'))))
kirc_cluster3_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster3_lls_hier.txt'))))
kirc_cluster2_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster4_lls_hier.txt'))))

kirc_cluster5_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster1_lls_iclust.txt'))))
kirc_cluster4_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster2_lls_iclust.txt'))))
kirc_cluster3_union_lls_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster3_lls_iclust.txt'))))
kirc_cluster2_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster4_lls_iclust.txt'))))

kirc_cluster5_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster1_lls_snf.txt'))))
kirc_cluster4_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster2_lls_snf.txt'))))
kirc_cluster3_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster3_lls_snf.txt'))))
kirc_cluster2_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
kirc_cluster5_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster1_lsa_hier.txt'))))
kirc_cluster4_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster2_lsa_hier.txt'))))
kirc_cluster3_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster3_lsa_hier.txt'))))
kirc_cluster2_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster4_lsa_hier.txt'))))

kirc_cluster5_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster1_lsa_iclust.txt'))))
kirc_cluster4_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster2_lsa_iclust.txt'))))
kirc_cluster3_union_lsa_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster3_lsa_iclust.txt'))))
kirc_cluster2_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster4_lsa_iclust.txt'))))

kirc_cluster5_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster1_lsa_snf.txt'))))
kirc_cluster4_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster2_lsa_snf.txt'))))
kirc_cluster3_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster3_lsa_snf.txt'))))
kirc_cluster2_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
kirc_cluster5_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster1_rand_hier.txt'))))
kirc_cluster4_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster2_rand_hier.txt'))))
kirc_cluster3_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster3_rand_hier.txt'))))
kirc_cluster2_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster4_rand_hier.txt'))))

kirc_cluster5_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster1_rand_iclust.txt'))))
kirc_cluster4_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster2_rand_iclust.txt'))))
kirc_cluster3_union_rand_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster3_rand_iclust.txt'))))
kirc_cluster2_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster4_rand_iclust.txt'))))

kirc_cluster5_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster1_rand_snf.txt'))))
kirc_cluster4_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster2_rand_snf.txt'))))
kirc_cluster3_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster3_rand_snf.txt'))))
kirc_cluster2_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/rand', '/cluster4_rand_snf.txt'))))


# similarity
kirc_cluster5_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster1_med_snf.txt'))))
kirc_cluster4_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster2_med_snf.txt'))))
kirc_cluster3_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster3_med_snf.txt'))))
kirc_cluster2_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster4_med_snf.txt'))))


kirc_cluster5_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster1_regress_snf.txt'))))
kirc_cluster4_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster2_regress_snf.txt'))))
kirc_cluster3_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster3_regress_snf.txt'))))
kirc_cluster2_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster4_regress_snf.txt'))))

kirc_cluster5_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster1_self_snf.txt'))))
kirc_cluster4_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster2_self_snf.txt'))))
kirc_cluster3_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster3_self_snf.txt'))))
kirc_cluster2_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/cluster4_self_snf.txt'))))


##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical############################################################################################

# Hierarchical 5 clusters
kirc_hier_stat_5 <- rbind(
  append('kirc_com_hier_5', summary(as.factor(kirc_com5_hier))),
  append('kirc_com_hier_knn_5', summary(as.factor(kirc_cluster5_com_knn_hier))),
  append('kirc_com_hier_lls_5', summary(as.factor(kirc_cluster5_com_lls_hier))),
  append('kirc_com_hier_lsa_5', summary(as.factor(kirc_cluster5_com_lsa_hier))),
  append('kirc_com_hier_rand_5', summary(as.factor(kirc_cluster5_com_rand_hier))),
  append('kirc_union_hier_knn_5', summary(as.factor(kirc_cluster5_union_knn_hier))),
  append('kirc_union_hier_lls_5', summary(as.factor(kirc_cluster5_union_lls_hier))),
  append('kirc_union_hier_lsa_5', summary(as.factor(kirc_cluster5_union_lsa_hier))),
  append('kirc_union_hier_rand_5', summary(as.factor(kirc_cluster5_union_rand_hier)))
)

# Hierarchical 4 clusters
kirc_hier_stat_4 <- rbind(
  append('kirc_com_hier_4', summary(as.factor(kirc_com4_hier))),
  append('kirc_com_hier_knn_4', summary(as.factor(kirc_cluster4_com_knn_hier))),
  append('kirc_com_hier_lls_4', summary(as.factor(kirc_cluster4_com_lls_hier))),
  append('kirc_com_hier_lsa_4', summary(as.factor(kirc_cluster4_com_lsa_hier))),
  append('kirc_com_hier_rand_4', summary(as.factor(kirc_cluster4_com_rand_hier))),
  append('kirc_union_hier_knn_4', summary(as.factor(kirc_cluster4_union_knn_hier))),
  append('kirc_union_hier_lls_4', summary(as.factor(kirc_cluster4_union_lls_hier))),
  append('kirc_union_hier_lsa_4', summary(as.factor(kirc_cluster4_union_lsa_hier))),
  append('kirc_union_hier_rand_4', summary(as.factor(kirc_cluster4_union_rand_hier)))
)

# Hierarchical 3 clusters
kirc_hier_stat_3 <- rbind(
  append('kirc_com_hier_3', summary(as.factor(kirc_com3_hier))),
  append('kirc_com_hier_knn_3', summary(as.factor(kirc_cluster3_com_knn_hier))),
  append('kirc_com_hier_lls_3', summary(as.factor(kirc_cluster3_com_lls_hier))),
  append('kirc_com_hier_lsa_3', summary(as.factor(kirc_cluster3_com_lsa_hier))),
  append('kirc_com_hier_rand_3', summary(as.factor(kirc_cluster3_com_rand_hier))),
  append('kirc_union_hier_knn_3', summary(as.factor(kirc_cluster3_union_knn_hier))),
  append('kirc_union_hier_lls_3', summary(as.factor(kirc_cluster3_union_lls_hier))),
  append('kirc_union_hier_lsa_3', summary(as.factor(kirc_cluster3_union_lsa_hier))),
  append('kirc_union_hier_rand_3', summary(as.factor(kirc_cluster3_union_rand_hier)))
)

# Hierarchical 2 clusters
kirc_hier_stat_2 <- rbind(
  append('kirc_com_hier_2', summary(as.factor(kirc_com2_hier))),
  append('kirc_com_hier_knn_2', summary(as.factor(kirc_cluster2_com_knn_hier))),
  append('kirc_com_hier_lls_2', summary(as.factor(kirc_cluster2_com_lls_hier))),
  append('kirc_com_hier_lsa_2', summary(as.factor(kirc_cluster2_com_lsa_hier))),
  append('kirc_com_hier_rand_2', summary(as.factor(kirc_cluster2_com_rand_hier))),
  append('kirc_union_hier_knn_2', summary(as.factor(kirc_cluster2_union_knn_hier))),
  append('kirc_union_hier_lls_2', summary(as.factor(kirc_cluster2_union_lls_hier))),
  append('kirc_union_hier_lsa_2', summary(as.factor(kirc_cluster2_union_lsa_hier))),
  append('kirc_union_hier_rand_2', summary(as.factor(kirc_cluster2_union_rand_hier)))
)


# icluster############################################################################################

# icluster 5 clusters
kirc_iclust_stat_5 <- rbind(
  append('kirc_com_iclust_5', summary(as.factor(kirc_com5_iclust))),
  append('kirc_com_iclust_knn_5', summary(as.factor(kirc_cluster5_com_knn_iclust))),
  append('kirc_com_iclust_lls_5', summary(as.factor(kirc_cluster5_com_lls_iclust))),
  append('kirc_com_iclust_lsa_5', summary(as.factor(kirc_cluster5_com_lsa_iclust))),
  append('kirc_com_iclust_rand_5', summary(as.factor(kirc_cluster5_com_rand_iclust))),
  append('kirc_union_iclust_knn_5', summary(as.factor(kirc_cluster5_union_knn_iclust))),
  append('kirc_union_iclust_lls_5', summary(as.factor(kirc_cluster5_union_lls_iclust))),
  append('kirc_union_iclust_lsa_5', summary(as.factor(kirc_cluster5_union_lsa_iclust))),
  append('kirc_union_iclust_rand_5', summary(as.factor(kirc_cluster5_union_rand_iclust)))
)

# icluster 4 clusters
kirc_iclust_stat_4 <- rbind(
  append('kirc_com_iclust_4', summary(as.factor(kirc_com4_iclust))),
  append('kirc_com_iclust_knn_4', summary(as.factor(kirc_cluster4_com_knn_iclust))),
  append('kirc_com_iclust_lls_4', summary(as.factor(kirc_cluster4_com_lls_iclust))),
  append('kirc_com_iclust_lsa_4', summary(as.factor(kirc_cluster4_com_lsa_iclust))),
  append('kirc_com_iclust_rand_4', summary(as.factor(kirc_cluster4_com_rand_iclust))),
  append('kirc_union_iclust_knn_4', summary(as.factor(kirc_cluster4_union_knn_iclust))),
  append('kirc_union_iclust_lls_4', summary(as.factor(kirc_cluster4_union_lls_iclust))),
  append('kirc_union_iclust_lsa_4', summary(as.factor(kirc_cluster4_union_lsa_iclust))),
  append('kirc_union_iclust_rand_4', summary(as.factor(kirc_cluster4_union_rand_iclust)))
)

# icluster 3 clusters
kirc_iclust_stat_3 <- rbind(
  append('kirc_com_iclust_3', summary(as.factor(kirc_com3_iclust))),
  append('kirc_com_iclust_knn_3', summary(as.factor(kirc_cluster3_com_knn_iclust))),
  append('kirc_com_iclust_lls_3', summary(as.factor(kirc_cluster3_com_lls_iclust))),
  append('kirc_com_iclust_lsa_3', summary(as.factor(kirc_cluster3_com_lsa_iclust))),
  append('kirc_com_iclust_rand_3', summary(as.factor(kirc_cluster3_com_rand_iclust))),
  append('kirc_union_iclust_knn_3', summary(as.factor(kirc_cluster3_union_knn_iclust))),
  append('kirc_union_iclust_lls_3', summary(as.factor(kirc_cluster3_union_lls_iclust))),
  append('kirc_union_iclust_lsa_3', summary(as.factor(kirc_cluster3_union_lsa_iclust))),
  append('kirc_union_iclust_rand_3', summary(as.factor(kirc_cluster3_union_rand_iclust)))
)

# icluster 2 clusters
kirc_iclust_stat_2 <- rbind(
  append('kirc_com_iclust_2', summary(as.factor(kirc_com2_iclust))),
  append('kirc_com_iclust_knn_2', summary(as.factor(kirc_cluster2_com_knn_iclust))),
  append('kirc_com_iclust_lls_2', summary(as.factor(kirc_cluster2_com_lls_iclust))),
  append('kirc_com_iclust_lsa_2', summary(as.factor(kirc_cluster2_com_lsa_iclust))),
  append('kirc_com_iclust_rand_2', summary(as.factor(kirc_cluster2_com_rand_iclust))),
  append('kirc_union_iclust_knn_2', summary(as.factor(kirc_cluster2_union_knn_iclust))),
  append('kirc_union_iclust_lls_2', summary(as.factor(kirc_cluster2_union_lls_iclust))),
  append('kirc_union_iclust_lsa_2', summary(as.factor(kirc_cluster2_union_lsa_iclust))),
  append('kirc_union_iclust_rand_2', summary(as.factor(kirc_cluster2_union_rand_iclust)))
)


# snf############################################################################################

# snf 5 clusters
kirc_snf_stat_5 <- rbind(
  
  append('kirc_com_snf_5', summary(as.factor(kirc_com5_snf))),
  append('kirc_com_snf_knn_5', summary(as.factor(kirc_cluster5_com_knn_snf))),
  append('kirc_com_snf_lls_5', summary(as.factor(kirc_cluster5_com_lls_snf))),
  append('kirc_com_snf_lsa_5', summary(as.factor(kirc_cluster5_com_lsa_snf))),
  append('kirc_com_snf_rand_5', summary(as.factor(kirc_cluster5_com_rand_snf))),
  append('kirc_com_snf_reg_5', summary(as.factor(kirc_cluster5_com_reg_snf))),
  append('kirc_com_snf_med_5', summary(as.factor(kirc_cluster5_com_med_snf))),
  append('kirc_com_snf_self_5', summary(as.factor(kirc_cluster5_com_self_snf))),
  append('kirc_union_snf_knn_5', summary(as.factor(kirc_cluster5_union_knn_snf))),
  append('kirc_union_snf_lls_5', summary(as.factor(kirc_cluster5_union_lls_snf))),
  append('kirc_union_snf_lsa_5', summary(as.factor(kirc_cluster5_union_lsa_snf))),
  append('kirc_union_snf_rand_5', summary(as.factor(kirc_cluster5_union_rand_snf))),
  append('kirc_union_snf_reg_5', summary(as.factor(kirc_cluster5_union_reg_snf))),
  append('kirc_union_snf_med_5', summary(as.factor(kirc_cluster5_union_med_snf))),
  append('kirc_union_snf_self_5', summary(as.factor(kirc_cluster5_union_self_snf)))
  
)

# snf 4 clusters
kirc_snf_stat_4 <- rbind(
  
  append('kirc_com_snf_4', summary(as.factor(kirc_com4_snf))),
  append('kirc_com_snf_knn_4', summary(as.factor(kirc_cluster4_com_knn_snf))),
  append('kirc_com_snf_lls_4', summary(as.factor(kirc_cluster4_com_lls_snf))),
  append('kirc_com_snf_lsa_4', summary(as.factor(kirc_cluster4_com_lsa_snf))),
  append('kirc_com_snf_rand_4', summary(as.factor(kirc_cluster4_com_rand_snf))),
  append('kirc_com_snf_reg_4', summary(as.factor(kirc_cluster4_com_reg_snf))),
  append('kirc_com_snf_med_4', summary(as.factor(kirc_cluster4_com_med_snf))),
  append('kirc_com_snf_self_4', summary(as.factor(kirc_cluster4_com_self_snf))),
  append('kirc_union_snf_knn_4', summary(as.factor(kirc_cluster4_union_knn_snf))),
  append('kirc_union_snf_lls_4', summary(as.factor(kirc_cluster4_union_lls_snf))),
  append('kirc_union_snf_lsa_4', summary(as.factor(kirc_cluster4_union_lsa_snf))),
  append('kirc_union_snf_rand_4', summary(as.factor(kirc_cluster4_union_rand_snf))),
  append('kirc_union_snf_reg_4', summary(as.factor(kirc_cluster4_union_reg_snf))),
  append('kirc_union_snf_med_4', summary(as.factor(kirc_cluster4_union_med_snf))),
  append('kirc_union_snf_self_4', summary(as.factor(kirc_cluster4_union_self_snf)))
  
)


# snf 3 clusters
kirc_snf_stat_3 <- rbind(
  
  append('kirc_com_snf_3', summary(as.factor(kirc_com3_snf))),
  append('kirc_com_snf_knn_3', summary(as.factor(kirc_cluster3_com_knn_snf))),
  append('kirc_com_snf_lls_3', summary(as.factor(kirc_cluster3_com_lls_snf))),
  append('kirc_com_snf_lsa_3', summary(as.factor(kirc_cluster3_com_lsa_snf))),
  append('kirc_com_snf_rand_3', summary(as.factor(kirc_cluster3_com_rand_snf))),
  append('kirc_com_snf_reg_3', summary(as.factor(kirc_cluster3_com_reg_snf))),
  append('kirc_com_snf_med_3', summary(as.factor(kirc_cluster3_com_med_snf))),
  append('kirc_com_snf_self_3', summary(as.factor(kirc_cluster3_com_self_snf))),
  append('kirc_union_snf_knn_3', summary(as.factor(kirc_cluster3_union_knn_snf))),
  append('kirc_union_snf_lls_3', summary(as.factor(kirc_cluster3_union_lls_snf))),
  append('kirc_union_snf_lsa_3', summary(as.factor(kirc_cluster3_union_lsa_snf))),
  append('kirc_union_snf_rand_3', summary(as.factor(kirc_cluster3_union_rand_snf))),
  append('kirc_union_snf_reg_3', summary(as.factor(kirc_cluster3_union_reg_snf))),
  append('kirc_union_snf_med_3', summary(as.factor(kirc_cluster3_union_med_snf))),
  append('kirc_union_snf_self_3', summary(as.factor(kirc_cluster3_union_self_snf)))
  
)

# snf 2 clusters
kirc_snf_stat_2 <- rbind(
  
  append('kirc_com_snf_2', summary(as.factor(kirc_com2_snf))),
  append('kirc_com_snf_knn_2', summary(as.factor(kirc_cluster2_com_knn_snf))),
  append('kirc_com_snf_lls_2', summary(as.factor(kirc_cluster2_com_lls_snf))),
  append('kirc_com_snf_lsa_2', summary(as.factor(kirc_cluster2_com_lsa_snf))),
  append('kirc_com_snf_rand_2', summary(as.factor(kirc_cluster2_com_rand_snf))),
  append('kirc_com_snf_reg_2', summary(as.factor(kirc_cluster2_com_reg_snf))),
  append('kirc_com_snf_med_2', summary(as.factor(kirc_cluster2_com_med_snf))),
  append('kirc_com_snf_self_2', summary(as.factor(kirc_cluster2_com_self_snf))),
  append('kirc_union_snf_knn_2', summary(as.factor(kirc_cluster2_union_knn_snf))),
  append('kirc_union_snf_lls_2', summary(as.factor(kirc_cluster2_union_lls_snf))),
  append('kirc_union_snf_lsa_2', summary(as.factor(kirc_cluster2_union_lsa_snf))),
  append('kirc_union_snf_rand_2', summary(as.factor(kirc_cluster2_union_rand_snf))),
  append('kirc_union_snf_reg_2', summary(as.factor(kirc_cluster2_union_reg_snf))),
  append('kirc_union_snf_med_2', summary(as.factor(kirc_cluster2_union_med_snf))),
  append('kirc_union_snf_self_2', summary(as.factor(kirc_cluster2_union_self_snf)))
  
)

#################################################################################################

##########################################################################################################
# lihc, compare clusters from complete, compete imputed, and union

# get lihc IDs
lihc_complete_ids <- lihc_ids[[1]]
lihc_union_ids <- lihc_ids[[2]]

# Get labels for each clustering type, from the intersection for lihc
lihc_com5_hier <- as.factor(t(read.table(paste0(completeFolder, '/3_1_1.txt'))))
lihc_com5_iclust <- as.factor(t(read.table(paste0(completeFolder, '/3_1_2.txt'))))
lihc_com5_snf <- as.factor(t(read.table(paste0(completeFolder, '/3_1_3.txt'))))

lihc_com4_hier <- as.factor(t(read.table(paste0(completeFolder, '/3_2_1.txt'))))
lihc_com4_iclust <- as.factor(t(read.table(paste0(completeFolder, '/3_2_2.txt'))))
lihc_com4_snf <- as.factor(t(read.table(paste0(completeFolder, '/3_2_3.txt'))))

lihc_com3_hier <- as.factor(t(read.table(paste0(completeFolder, '/3_3_1.txt'))))
lihc_com3_iclust <- as.factor(t(read.table(paste0(completeFolder, '/3_3_2.txt'))))
lihc_com3_snf <- as.factor(t(read.table(paste0(completeFolder, '/3_3_3.txt'))))

lihc_com2_hier <- as.factor(t(read.table(paste0(completeFolder, '/3_4_1.txt'))))
lihc_com2_iclust <- as.factor(t(read.table(paste0(completeFolder, '/3_4_2.txt'))))
lihc_com2_snf <- as.factor(t(read.table(paste0(completeFolder, '/3_4_3.txt'))))
############################################################################################3
# get labels for lihc intersection imputed
# knn cluster 1 - 4
lihc_cluster5_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster1_knn_hier.txt'))))
lihc_cluster4_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster2_knn_hier.txt'))))
lihc_cluster3_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster3_knn_hier.txt'))))
lihc_cluster2_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster4_knn_hier.txt'))))

lihc_cluster5_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster1_knn_iclust.txt'))))
lihc_cluster4_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster2_knn_iclust.txt'))))
lihc_cluster3_com_knn_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster3_knn_iclust.txt'))))
lihc_cluster2_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster4_knn_iclust.txt'))))

lihc_cluster5_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster1_knn_snf.txt'))))
lihc_cluster4_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster2_knn_snf.txt'))))
lihc_cluster3_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster3_knn_snf.txt'))))
lihc_cluster2_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
lihc_cluster5_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster1_lls_hier.txt'))))
lihc_cluster4_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster2_lls_hier.txt'))))
lihc_cluster3_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster3_lls_hier.txt'))))
lihc_cluster2_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster4_lls_hier.txt'))))

lihc_cluster5_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster1_lls_iclust.txt'))))
lihc_cluster4_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster2_lls_iclust.txt'))))
lihc_cluster3_com_lls_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster3_lls_iclust.txt'))))
lihc_cluster2_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster4_lls_iclust.txt'))))

lihc_cluster5_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster1_lls_snf.txt'))))
lihc_cluster4_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster2_lls_snf.txt'))))
lihc_cluster3_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster3_lls_snf.txt'))))
lihc_cluster2_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
lihc_cluster5_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster1_lsa_hier.txt'))))
lihc_cluster4_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster2_lsa_hier.txt'))))
lihc_cluster3_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster3_lsa_hier.txt'))))
lihc_cluster2_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster4_lsa_hier.txt'))))

lihc_cluster5_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster1_lsa_iclust.txt'))))
lihc_cluster4_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster2_lsa_iclust.txt'))))
lihc_cluster3_com_lsa_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster3_lsa_iclust.txt'))))
lihc_cluster2_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster4_lsa_iclust.txt'))))

lihc_cluster5_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster1_lsa_snf.txt'))))
lihc_cluster4_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster2_lsa_snf.txt'))))
lihc_cluster3_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster3_lsa_snf.txt'))))
lihc_cluster2_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
lihc_cluster5_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster1_rand_hier.txt'))))
lihc_cluster4_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster2_rand_hier.txt'))))
lihc_cluster3_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster3_rand_hier.txt'))))
lihc_cluster2_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster4_rand_hier.txt'))))

lihc_cluster5_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster1_rand_iclust.txt'))))
lihc_cluster4_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster2_rand_iclust.txt'))))
lihc_cluster3_com_rand_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster3_rand_iclust.txt'))))
lihc_cluster2_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster4_rand_iclust.txt'))))

lihc_cluster5_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster1_rand_snf.txt'))))
lihc_cluster4_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster2_rand_snf.txt'))))
lihc_cluster3_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster3_rand_snf.txt'))))
lihc_cluster2_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lihc/rand', '/cluster4_rand_snf.txt'))))


# similarity
lihc_cluster5_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster1_med_snf.txt'))))
lihc_cluster4_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster2_med_snf.txt'))))
lihc_cluster3_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster3_med_snf.txt'))))
lihc_cluster2_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster4_med_snf.txt'))))


lihc_cluster5_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster1_regress_snf.txt'))))
lihc_cluster4_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster2_regress_snf.txt'))))
lihc_cluster3_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster3_regress_snf.txt'))))
lihc_cluster2_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster4_regress_snf.txt'))))

lihc_cluster5_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster1_self_snf.txt'))))
lihc_cluster4_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster2_self_snf.txt'))))
lihc_cluster3_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster3_self_snf.txt'))))
lihc_cluster2_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/cluster4_self_snf.txt'))))

#########################################################################################################
# Get labels for union

# knn cluster 1 - 4
lihc_cluster5_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster1_knn_hier.txt'))))
lihc_cluster4_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster2_knn_hier.txt'))))
lihc_cluster3_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster3_knn_hier.txt'))))
lihc_cluster2_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster4_knn_hier.txt'))))

lihc_cluster5_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster1_knn_iclust.txt'))))
lihc_cluster4_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster2_knn_iclust.txt'))))
lihc_cluster3_union_knn_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster3_knn_iclust.txt'))))
lihc_cluster2_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster4_knn_iclust.txt'))))

lihc_cluster5_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster1_knn_snf.txt'))))
lihc_cluster4_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster2_knn_snf.txt'))))
lihc_cluster3_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster3_knn_snf.txt'))))
lihc_cluster2_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
lihc_cluster5_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster1_lls_hier.txt'))))
lihc_cluster4_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster2_lls_hier.txt'))))
lihc_cluster3_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster3_lls_hier.txt'))))
lihc_cluster2_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster4_lls_hier.txt'))))

lihc_cluster5_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster1_lls_iclust.txt'))))
lihc_cluster4_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster2_lls_iclust.txt'))))
lihc_cluster3_union_lls_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster3_lls_iclust.txt'))))
lihc_cluster2_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster4_lls_iclust.txt'))))

lihc_cluster5_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster1_lls_snf.txt'))))
lihc_cluster4_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster2_lls_snf.txt'))))
lihc_cluster3_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster3_lls_snf.txt'))))
lihc_cluster2_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
lihc_cluster5_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster1_lsa_hier.txt'))))
lihc_cluster4_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster2_lsa_hier.txt'))))
lihc_cluster3_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster3_lsa_hier.txt'))))
lihc_cluster2_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster4_lsa_hier.txt'))))

lihc_cluster5_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster1_lsa_iclust.txt'))))
lihc_cluster4_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster2_lsa_iclust.txt'))))
lihc_cluster3_union_lsa_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster3_lsa_iclust.txt'))))
lihc_cluster2_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster4_lsa_iclust.txt'))))

lihc_cluster5_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster1_lsa_snf.txt'))))
lihc_cluster4_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster2_lsa_snf.txt'))))
lihc_cluster3_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster3_lsa_snf.txt'))))
lihc_cluster2_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
lihc_cluster5_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster1_rand_hier.txt'))))
lihc_cluster4_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster2_rand_hier.txt'))))
lihc_cluster3_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster3_rand_hier.txt'))))
lihc_cluster2_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster4_rand_hier.txt'))))

lihc_cluster5_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster1_rand_iclust.txt'))))
lihc_cluster4_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster2_rand_iclust.txt'))))
lihc_cluster3_union_rand_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster3_rand_iclust.txt'))))
lihc_cluster2_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster4_rand_iclust.txt'))))

lihc_cluster5_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster1_rand_snf.txt'))))
lihc_cluster4_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster2_rand_snf.txt'))))
lihc_cluster3_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster3_rand_snf.txt'))))
lihc_cluster2_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/rand', '/cluster4_rand_snf.txt'))))


# similarity
lihc_cluster5_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster1_med_snf.txt'))))
lihc_cluster4_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster2_med_snf.txt'))))
lihc_cluster3_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster3_med_snf.txt'))))
lihc_cluster2_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster4_med_snf.txt'))))


lihc_cluster5_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster1_regress_snf.txt'))))
lihc_cluster4_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster2_regress_snf.txt'))))
lihc_cluster3_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster3_regress_snf.txt'))))
lihc_cluster2_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster4_regress_snf.txt'))))

lihc_cluster5_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster1_self_snf.txt'))))
lihc_cluster4_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster2_self_snf.txt'))))
lihc_cluster3_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster3_self_snf.txt'))))
lihc_cluster2_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/cluster4_self_snf.txt'))))


##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical############################################################################################

# Hierarchical 5 clusters
lihc_hier_stat_5 <- rbind(
  append('lihc_com_hier_5', summary(as.factor(lihc_com5_hier))),
  append('lihc_com_hier_knn_5', summary(as.factor(lihc_cluster5_com_knn_hier))),
  append('lihc_com_hier_lls_5', summary(as.factor(lihc_cluster5_com_lls_hier))),
  append('lihc_com_hier_lsa_5', summary(as.factor(lihc_cluster5_com_lsa_hier))),
  append('lihc_com_hier_rand_5', summary(as.factor(lihc_cluster5_com_rand_hier))),
  append('lihc_union_hier_knn_5', summary(as.factor(lihc_cluster5_union_knn_hier))),
  append('lihc_union_hier_lls_5', summary(as.factor(lihc_cluster5_union_lls_hier))),
  append('lihc_union_hier_lsa_5', summary(as.factor(lihc_cluster5_union_lsa_hier))),
  append('lihc_union_hier_rand_5', summary(as.factor(lihc_cluster5_union_rand_hier)))
)

# Hierarchical 4 clusters
lihc_hier_stat_4 <- rbind(
  append('lihc_com_hier_4', summary(as.factor(lihc_com4_hier))),
  append('lihc_com_hier_knn_4', summary(as.factor(lihc_cluster4_com_knn_hier))),
  append('lihc_com_hier_lls_4', summary(as.factor(lihc_cluster4_com_lls_hier))),
  append('lihc_com_hier_lsa_4', summary(as.factor(lihc_cluster4_com_lsa_hier))),
  append('lihc_com_hier_rand_4', summary(as.factor(lihc_cluster4_com_rand_hier))),
  append('lihc_union_hier_knn_4', summary(as.factor(lihc_cluster4_union_knn_hier))),
  append('lihc_union_hier_lls_4', summary(as.factor(lihc_cluster4_union_lls_hier))),
  append('lihc_union_hier_lsa_4', summary(as.factor(lihc_cluster4_union_lsa_hier))),
  append('lihc_union_hier_rand_4', summary(as.factor(lihc_cluster4_union_rand_hier)))
)

# Hierarchical 3 clusters
lihc_hier_stat_3 <- rbind(
  append('lihc_com_hier_3', summary(as.factor(lihc_com3_hier))),
  append('lihc_com_hier_knn_3', summary(as.factor(lihc_cluster3_com_knn_hier))),
  append('lihc_com_hier_lls_3', summary(as.factor(lihc_cluster3_com_lls_hier))),
  append('lihc_com_hier_lsa_3', summary(as.factor(lihc_cluster3_com_lsa_hier))),
  append('lihc_com_hier_rand_3', summary(as.factor(lihc_cluster3_com_rand_hier))),
  append('lihc_union_hier_knn_3', summary(as.factor(lihc_cluster3_union_knn_hier))),
  append('lihc_union_hier_lls_3', summary(as.factor(lihc_cluster3_union_lls_hier))),
  append('lihc_union_hier_lsa_3', summary(as.factor(lihc_cluster3_union_lsa_hier))),
  append('lihc_union_hier_rand_3', summary(as.factor(lihc_cluster3_union_rand_hier)))
)

# Hierarchical 2 clusters
lihc_hier_stat_2 <- rbind(
  append('lihc_com_hier_2', summary(as.factor(lihc_com2_hier))),
  append('lihc_com_hier_knn_2', summary(as.factor(lihc_cluster2_com_knn_hier))),
  append('lihc_com_hier_lls_2', summary(as.factor(lihc_cluster2_com_lls_hier))),
  append('lihc_com_hier_lsa_2', summary(as.factor(lihc_cluster2_com_lsa_hier))),
  append('lihc_com_hier_rand_2', summary(as.factor(lihc_cluster2_com_rand_hier))),
  append('lihc_union_hier_knn_2', summary(as.factor(lihc_cluster2_union_knn_hier))),
  append('lihc_union_hier_lls_2', summary(as.factor(lihc_cluster2_union_lls_hier))),
  append('lihc_union_hier_lsa_2', summary(as.factor(lihc_cluster2_union_lsa_hier))),
  append('lihc_union_hier_rand_2', summary(as.factor(lihc_cluster2_union_rand_hier)))
)


# icluster############################################################################################

# icluster 5 clusters
lihc_iclust_stat_5 <- rbind(
  append('lihc_com_iclust_5', summary(as.factor(lihc_com5_iclust))),
  append('lihc_com_iclust_knn_5', summary(as.factor(lihc_cluster5_com_knn_iclust))),
  append('lihc_com_iclust_lls_5', summary(as.factor(lihc_cluster5_com_lls_iclust))),
  append('lihc_com_iclust_lsa_5', summary(as.factor(lihc_cluster5_com_lsa_iclust))),
  append('lihc_com_iclust_rand_5', summary(as.factor(lihc_cluster5_com_rand_iclust))),
  append('lihc_union_iclust_knn_5', summary(as.factor(lihc_cluster5_union_knn_iclust))),
  append('lihc_union_iclust_lls_5', summary(as.factor(lihc_cluster5_union_lls_iclust))),
  append('lihc_union_iclust_lsa_5', summary(as.factor(lihc_cluster5_union_lsa_iclust))),
  append('lihc_union_iclust_rand_5', summary(as.factor(lihc_cluster5_union_rand_iclust)))
)

# icluster 4 clusters
lihc_iclust_stat_4 <- rbind(
  append('lihc_com_iclust_4', summary(as.factor(lihc_com4_iclust))),
  append('lihc_com_iclust_knn_4', summary(as.factor(lihc_cluster4_com_knn_iclust))),
  append('lihc_com_iclust_lls_4', summary(as.factor(lihc_cluster4_com_lls_iclust))),
  append('lihc_com_iclust_lsa_4', summary(as.factor(lihc_cluster4_com_lsa_iclust))),
  append('lihc_com_iclust_rand_4', summary(as.factor(lihc_cluster4_com_rand_iclust))),
  append('lihc_union_iclust_knn_4', summary(as.factor(lihc_cluster4_union_knn_iclust))),
  append('lihc_union_iclust_lls_4', summary(as.factor(lihc_cluster4_union_lls_iclust))),
  append('lihc_union_iclust_lsa_4', summary(as.factor(lihc_cluster4_union_lsa_iclust))),
  append('lihc_union_iclust_rand_4', summary(as.factor(lihc_cluster4_union_rand_iclust)))
)

# icluster 3 clusters
lihc_iclust_stat_3 <- rbind(
  append('lihc_com_iclust_3', summary(as.factor(lihc_com3_iclust))),
  append('lihc_com_iclust_knn_3', summary(as.factor(lihc_cluster3_com_knn_iclust))),
  append('lihc_com_iclust_lls_3', summary(as.factor(lihc_cluster3_com_lls_iclust))),
  append('lihc_com_iclust_lsa_3', summary(as.factor(lihc_cluster3_com_lsa_iclust))),
  append('lihc_com_iclust_rand_3', summary(as.factor(lihc_cluster3_com_rand_iclust))),
  append('lihc_union_iclust_knn_3', summary(as.factor(lihc_cluster3_union_knn_iclust))),
  append('lihc_union_iclust_lls_3', summary(as.factor(lihc_cluster3_union_lls_iclust))),
  append('lihc_union_iclust_lsa_3', summary(as.factor(lihc_cluster3_union_lsa_iclust))),
  append('lihc_union_iclust_rand_3', summary(as.factor(lihc_cluster3_union_rand_iclust)))
)

# icluster 2 clusters
lihc_iclust_stat_2 <- rbind(
  append('lihc_com_iclust_2', summary(as.factor(lihc_com2_iclust))),
  append('lihc_com_iclust_knn_2', summary(as.factor(lihc_cluster2_com_knn_iclust))),
  append('lihc_com_iclust_lls_2', summary(as.factor(lihc_cluster2_com_lls_iclust))),
  append('lihc_com_iclust_lsa_2', summary(as.factor(lihc_cluster2_com_lsa_iclust))),
  append('lihc_com_iclust_rand_2', summary(as.factor(lihc_cluster2_com_rand_iclust))),
  append('lihc_union_iclust_knn_2', summary(as.factor(lihc_cluster2_union_knn_iclust))),
  append('lihc_union_iclust_lls_2', summary(as.factor(lihc_cluster2_union_lls_iclust))),
  append('lihc_union_iclust_lsa_2', summary(as.factor(lihc_cluster2_union_lsa_iclust))),
  append('lihc_union_iclust_rand_2', summary(as.factor(lihc_cluster2_union_rand_iclust)))
)


# snf############################################################################################

# snf 5 clusters
lihc_snf_stat_5 <- rbind(
  
  append('lihc_com_snf_5', summary(as.factor(lihc_com5_snf))),
  append('lihc_com_snf_knn_5', summary(as.factor(lihc_cluster5_com_knn_snf))),
  append('lihc_com_snf_lls_5', summary(as.factor(lihc_cluster5_com_lls_snf))),
  append('lihc_com_snf_lsa_5', summary(as.factor(lihc_cluster5_com_lsa_snf))),
  append('lihc_com_snf_rand_5', summary(as.factor(lihc_cluster5_com_rand_snf))),
  append('lihc_com_snf_reg_5', summary(as.factor(lihc_cluster5_com_reg_snf))),
  append('lihc_com_snf_med_5', summary(as.factor(lihc_cluster5_com_med_snf))),
  append('lihc_com_snf_self_5', summary(as.factor(lihc_cluster5_com_self_snf))),
  append('lihc_union_snf_knn_5', summary(as.factor(lihc_cluster5_union_knn_snf))),
  append('lihc_union_snf_lls_5', summary(as.factor(lihc_cluster5_union_lls_snf))),
  append('lihc_union_snf_lsa_5', summary(as.factor(lihc_cluster5_union_lsa_snf))),
  append('lihc_union_snf_rand_5', summary(as.factor(lihc_cluster5_union_rand_snf))),
  append('lihc_union_snf_reg_5', summary(as.factor(lihc_cluster5_union_reg_snf))),
  append('lihc_union_snf_med_5', summary(as.factor(lihc_cluster5_union_med_snf))),
  append('lihc_union_snf_self_5', summary(as.factor(lihc_cluster5_union_self_snf)))
  
)

# snf 4 clusters
lihc_snf_stat_4 <- rbind(
  
  append('lihc_com_snf_4', summary(as.factor(lihc_com4_snf))),
  append('lihc_com_snf_knn_4', summary(as.factor(lihc_cluster4_com_knn_snf))),
  append('lihc_com_snf_lls_4', summary(as.factor(lihc_cluster4_com_lls_snf))),
  append('lihc_com_snf_lsa_4', summary(as.factor(lihc_cluster4_com_lsa_snf))),
  append('lihc_com_snf_rand_4', summary(as.factor(lihc_cluster4_com_rand_snf))),
  append('lihc_com_snf_reg_4', summary(as.factor(lihc_cluster4_com_reg_snf))),
  append('lihc_com_snf_med_4', summary(as.factor(lihc_cluster4_com_med_snf))),
  append('lihc_com_snf_self_4', summary(as.factor(lihc_cluster4_com_self_snf))),
  append('lihc_union_snf_knn_4', summary(as.factor(lihc_cluster4_union_knn_snf))),
  append('lihc_union_snf_lls_4', summary(as.factor(lihc_cluster4_union_lls_snf))),
  append('lihc_union_snf_lsa_4', summary(as.factor(lihc_cluster4_union_lsa_snf))),
  append('lihc_union_snf_rand_4', summary(as.factor(lihc_cluster4_union_rand_snf))),
  append('lihc_union_snf_reg_4', summary(as.factor(lihc_cluster4_union_reg_snf))),
  append('lihc_union_snf_med_4', summary(as.factor(lihc_cluster4_union_med_snf))),
  append('lihc_union_snf_self_4', summary(as.factor(lihc_cluster4_union_self_snf)))
  
)


# snf 3 clusters
lihc_snf_stat_3 <- rbind(
  
  append('lihc_com_snf_3', summary(as.factor(lihc_com3_snf))),
  append('lihc_com_snf_knn_3', summary(as.factor(lihc_cluster3_com_knn_snf))),
  append('lihc_com_snf_lls_3', summary(as.factor(lihc_cluster3_com_lls_snf))),
  append('lihc_com_snf_lsa_3', summary(as.factor(lihc_cluster3_com_lsa_snf))),
  append('lihc_com_snf_rand_3', summary(as.factor(lihc_cluster3_com_rand_snf))),
  append('lihc_com_snf_reg_3', summary(as.factor(lihc_cluster3_com_reg_snf))),
  append('lihc_com_snf_med_3', summary(as.factor(lihc_cluster3_com_med_snf))),
  append('lihc_com_snf_self_3', summary(as.factor(lihc_cluster3_com_self_snf))),
  append('lihc_union_snf_knn_3', summary(as.factor(lihc_cluster3_union_knn_snf))),
  append('lihc_union_snf_lls_3', summary(as.factor(lihc_cluster3_union_lls_snf))),
  append('lihc_union_snf_lsa_3', summary(as.factor(lihc_cluster3_union_lsa_snf))),
  append('lihc_union_snf_rand_3', summary(as.factor(lihc_cluster3_union_rand_snf))),
  append('lihc_union_snf_reg_3', summary(as.factor(lihc_cluster3_union_reg_snf))),
  append('lihc_union_snf_med_3', summary(as.factor(lihc_cluster3_union_med_snf))),
  append('lihc_union_snf_self_3', summary(as.factor(lihc_cluster3_union_self_snf)))
  
)

# snf 2 clusters
lihc_snf_stat_2 <- rbind(
  
  append('lihc_com_snf_2', summary(as.factor(lihc_com2_snf))),
  append('lihc_com_snf_knn_2', summary(as.factor(lihc_cluster2_com_knn_snf))),
  append('lihc_com_snf_lls_2', summary(as.factor(lihc_cluster2_com_lls_snf))),
  append('lihc_com_snf_lsa_2', summary(as.factor(lihc_cluster2_com_lsa_snf))),
  append('lihc_com_snf_rand_2', summary(as.factor(lihc_cluster2_com_rand_snf))),
  append('lihc_com_snf_reg_2', summary(as.factor(lihc_cluster2_com_reg_snf))),
  append('lihc_com_snf_med_2', summary(as.factor(lihc_cluster2_com_med_snf))),
  append('lihc_com_snf_self_2', summary(as.factor(lihc_cluster2_com_self_snf))),
  append('lihc_union_snf_knn_2', summary(as.factor(lihc_cluster2_union_knn_snf))),
  append('lihc_union_snf_lls_2', summary(as.factor(lihc_cluster2_union_lls_snf))),
  append('lihc_union_snf_lsa_2', summary(as.factor(lihc_cluster2_union_lsa_snf))),
  append('lihc_union_snf_rand_2', summary(as.factor(lihc_cluster2_union_rand_snf))),
  append('lihc_union_snf_reg_2', summary(as.factor(lihc_cluster2_union_reg_snf))),
  append('lihc_union_snf_med_2', summary(as.factor(lihc_cluster2_union_med_snf))),
  append('lihc_union_snf_self_2', summary(as.factor(lihc_cluster2_union_self_snf)))
  
)

#################################################################################################

##########################################################################################################
# luad, compare clusters from complete, compete imputed, and union

# get luad IDs
luad_complete_ids <- luad_ids[[1]]
luad_union_ids <- luad_ids[[2]]

# Get labels for each clustering type, from the intersection for luad
luad_com5_hier <- as.factor(t(read.table(paste0(completeFolder, '/4_1_1.txt'))))
luad_com5_iclust <- as.factor(t(read.table(paste0(completeFolder, '/4_1_2.txt'))))
luad_com5_snf <- as.factor(t(read.table(paste0(completeFolder, '/4_1_3.txt'))))

luad_com4_hier <- as.factor(t(read.table(paste0(completeFolder, '/4_2_1.txt'))))
luad_com4_iclust <- as.factor(t(read.table(paste0(completeFolder, '/4_2_2.txt'))))
luad_com4_snf <- as.factor(t(read.table(paste0(completeFolder, '/4_2_3.txt'))))

luad_com3_hier <- as.factor(t(read.table(paste0(completeFolder, '/4_3_1.txt'))))
luad_com3_iclust <- as.factor(t(read.table(paste0(completeFolder, '/4_3_2.txt'))))
luad_com3_snf <- as.factor(t(read.table(paste0(completeFolder, '/4_3_3.txt'))))

luad_com2_hier <- as.factor(t(read.table(paste0(completeFolder, '/4_4_1.txt'))))
luad_com2_iclust <- as.factor(t(read.table(paste0(completeFolder, '/4_4_2.txt'))))
luad_com2_snf <- as.factor(t(read.table(paste0(completeFolder, '/4_4_3.txt'))))
############################################################################################3
# get labels for luad intersection imputed
# knn cluster 1 - 4
luad_cluster5_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster1_knn_hier.txt'))))
luad_cluster4_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster2_knn_hier.txt'))))
luad_cluster3_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster3_knn_hier.txt'))))
luad_cluster2_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster4_knn_hier.txt'))))

luad_cluster5_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster1_knn_iclust.txt'))))
luad_cluster4_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster2_knn_iclust.txt'))))
luad_cluster3_com_knn_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster3_knn_iclust.txt'))))
luad_cluster2_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster4_knn_iclust.txt'))))

luad_cluster5_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster1_knn_snf.txt'))))
luad_cluster4_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster2_knn_snf.txt'))))
luad_cluster3_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster3_knn_snf.txt'))))
luad_cluster2_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
luad_cluster5_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster1_lls_hier.txt'))))
luad_cluster4_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster2_lls_hier.txt'))))
luad_cluster3_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster3_lls_hier.txt'))))
luad_cluster2_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster4_lls_hier.txt'))))

luad_cluster5_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster1_lls_iclust.txt'))))
luad_cluster4_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster2_lls_iclust.txt'))))
luad_cluster3_com_lls_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster3_lls_iclust.txt'))))
luad_cluster2_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster4_lls_iclust.txt'))))

luad_cluster5_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster1_lls_snf.txt'))))
luad_cluster4_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster2_lls_snf.txt'))))
luad_cluster3_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster3_lls_snf.txt'))))
luad_cluster2_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
luad_cluster5_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster1_lsa_hier.txt'))))
luad_cluster4_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster2_lsa_hier.txt'))))
luad_cluster3_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster3_lsa_hier.txt'))))
luad_cluster2_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster4_lsa_hier.txt'))))

luad_cluster5_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster1_lsa_iclust.txt'))))
luad_cluster4_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster2_lsa_iclust.txt'))))
luad_cluster3_com_lsa_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster3_lsa_iclust.txt'))))
luad_cluster2_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster4_lsa_iclust.txt'))))

luad_cluster5_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster1_lsa_snf.txt'))))
luad_cluster4_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster2_lsa_snf.txt'))))
luad_cluster3_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster3_lsa_snf.txt'))))
luad_cluster2_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
luad_cluster5_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster1_rand_hier.txt'))))
luad_cluster4_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster2_rand_hier.txt'))))
luad_cluster3_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster3_rand_hier.txt'))))
luad_cluster2_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster4_rand_hier.txt'))))

luad_cluster5_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster1_rand_iclust.txt'))))
luad_cluster4_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster2_rand_iclust.txt'))))
luad_cluster3_com_rand_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster3_rand_iclust.txt'))))
luad_cluster2_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster4_rand_iclust.txt'))))

luad_cluster5_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster1_rand_snf.txt'))))
luad_cluster4_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster2_rand_snf.txt'))))
luad_cluster3_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster3_rand_snf.txt'))))
luad_cluster2_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/luad/rand', '/cluster4_rand_snf.txt'))))


# similarity
luad_cluster5_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster1_med_snf.txt'))))
luad_cluster4_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster2_med_snf.txt'))))
luad_cluster3_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster3_med_snf.txt'))))
luad_cluster2_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster4_med_snf.txt'))))


luad_cluster5_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster1_regress_snf.txt'))))
luad_cluster4_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster2_regress_snf.txt'))))
luad_cluster3_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster3_regress_snf.txt'))))
luad_cluster2_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster4_regress_snf.txt'))))

luad_cluster5_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster1_self_snf.txt'))))
luad_cluster4_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster2_self_snf.txt'))))
luad_cluster3_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster3_self_snf.txt'))))
luad_cluster2_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/cluster4_self_snf.txt'))))

#########################################################################################################
# Get labels for union

# knn cluster 1 - 4
luad_cluster5_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster1_knn_hier.txt'))))
luad_cluster4_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster2_knn_hier.txt'))))
luad_cluster3_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster3_knn_hier.txt'))))
luad_cluster2_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster4_knn_hier.txt'))))

luad_cluster5_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster1_knn_iclust.txt'))))
luad_cluster4_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster2_knn_iclust.txt'))))
luad_cluster3_union_knn_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster3_knn_iclust.txt'))))
luad_cluster2_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster4_knn_iclust.txt'))))

luad_cluster5_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster1_knn_snf.txt'))))
luad_cluster4_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster2_knn_snf.txt'))))
luad_cluster3_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster3_knn_snf.txt'))))
luad_cluster2_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
luad_cluster5_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster1_lls_hier.txt'))))
luad_cluster4_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster2_lls_hier.txt'))))
luad_cluster3_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster3_lls_hier.txt'))))
luad_cluster2_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster4_lls_hier.txt'))))

luad_cluster5_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster1_lls_iclust.txt'))))
luad_cluster4_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster2_lls_iclust.txt'))))
luad_cluster3_union_lls_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster3_lls_iclust.txt'))))
luad_cluster2_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster4_lls_iclust.txt'))))

luad_cluster5_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster1_lls_snf.txt'))))
luad_cluster4_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster2_lls_snf.txt'))))
luad_cluster3_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster3_lls_snf.txt'))))
luad_cluster2_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
luad_cluster5_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster1_lsa_hier.txt'))))
luad_cluster4_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster2_lsa_hier.txt'))))
luad_cluster3_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster3_lsa_hier.txt'))))
luad_cluster2_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster4_lsa_hier.txt'))))

luad_cluster5_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster1_lsa_iclust.txt'))))
luad_cluster4_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster2_lsa_iclust.txt'))))
luad_cluster3_union_lsa_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster3_lsa_iclust.txt'))))
luad_cluster2_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster4_lsa_iclust.txt'))))

luad_cluster5_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster1_lsa_snf.txt'))))
luad_cluster4_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster2_lsa_snf.txt'))))
luad_cluster3_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster3_lsa_snf.txt'))))
luad_cluster2_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
luad_cluster5_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster1_rand_hier.txt'))))
luad_cluster4_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster2_rand_hier.txt'))))
luad_cluster3_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster3_rand_hier.txt'))))
luad_cluster2_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster4_rand_hier.txt'))))

luad_cluster5_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster1_rand_iclust.txt'))))
luad_cluster4_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster2_rand_iclust.txt'))))
luad_cluster3_union_rand_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster3_rand_iclust.txt'))))
luad_cluster2_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster4_rand_iclust.txt'))))

luad_cluster5_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster1_rand_snf.txt'))))
luad_cluster4_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster2_rand_snf.txt'))))
luad_cluster3_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster3_rand_snf.txt'))))
luad_cluster2_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/rand', '/cluster4_rand_snf.txt'))))


# similarity
luad_cluster5_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster1_med_snf.txt'))))
luad_cluster4_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster2_med_snf.txt'))))
luad_cluster3_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster3_med_snf.txt'))))
luad_cluster2_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster4_med_snf.txt'))))


luad_cluster5_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster1_regress_snf.txt'))))
luad_cluster4_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster2_regress_snf.txt'))))
luad_cluster3_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster3_regress_snf.txt'))))
luad_cluster2_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster4_regress_snf.txt'))))

luad_cluster5_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster1_self_snf.txt'))))
luad_cluster4_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster2_self_snf.txt'))))
luad_cluster3_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster3_self_snf.txt'))))
luad_cluster2_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/cluster4_self_snf.txt'))))


##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical############################################################################################

# Hierarchical 5 clusters
luad_hier_stat_5 <- rbind(
  append('luad_com_hier_5', summary(as.factor(luad_com5_hier))),
  append('luad_com_hier_knn_5', summary(as.factor(luad_cluster5_com_knn_hier))),
  append('luad_com_hier_lls_5', summary(as.factor(luad_cluster5_com_lls_hier))),
  append('luad_com_hier_lsa_5', summary(as.factor(luad_cluster5_com_lsa_hier))),
  append('luad_com_hier_rand_5', summary(as.factor(luad_cluster5_com_rand_hier))),
  append('luad_union_hier_knn_5', summary(as.factor(luad_cluster5_union_knn_hier))),
  append('luad_union_hier_lls_5', summary(as.factor(luad_cluster5_union_lls_hier))),
  append('luad_union_hier_lsa_5', summary(as.factor(luad_cluster5_union_lsa_hier))),
  append('luad_union_hier_rand_5', summary(as.factor(luad_cluster5_union_rand_hier)))
)

# Hierarchical 4 clusters
luad_hier_stat_4 <- rbind(
  append('luad_com_hier_4', summary(as.factor(luad_com4_hier))),
  append('luad_com_hier_knn_4', summary(as.factor(luad_cluster4_com_knn_hier))),
  append('luad_com_hier_lls_4', summary(as.factor(luad_cluster4_com_lls_hier))),
  append('luad_com_hier_lsa_4', summary(as.factor(luad_cluster4_com_lsa_hier))),
  append('luad_com_hier_rand_4', summary(as.factor(luad_cluster4_com_rand_hier))),
  append('luad_union_hier_knn_4', summary(as.factor(luad_cluster4_union_knn_hier))),
  append('luad_union_hier_lls_4', summary(as.factor(luad_cluster4_union_lls_hier))),
  append('luad_union_hier_lsa_4', summary(as.factor(luad_cluster4_union_lsa_hier))),
  append('luad_union_hier_rand_4', summary(as.factor(luad_cluster4_union_rand_hier)))
)

# Hierarchical 3 clusters
luad_hier_stat_3 <- rbind(
  append('luad_com_hier_3', summary(as.factor(luad_com3_hier))),
  append('luad_com_hier_knn_3', summary(as.factor(luad_cluster3_com_knn_hier))),
  append('luad_com_hier_lls_3', summary(as.factor(luad_cluster3_com_lls_hier))),
  append('luad_com_hier_lsa_3', summary(as.factor(luad_cluster3_com_lsa_hier))),
  append('luad_com_hier_rand_3', summary(as.factor(luad_cluster3_com_rand_hier))),
  append('luad_union_hier_knn_3', summary(as.factor(luad_cluster3_union_knn_hier))),
  append('luad_union_hier_lls_3', summary(as.factor(luad_cluster3_union_lls_hier))),
  append('luad_union_hier_lsa_3', summary(as.factor(luad_cluster3_union_lsa_hier))),
  append('luad_union_hier_rand_3', summary(as.factor(luad_cluster3_union_rand_hier)))
)

# Hierarchical 2 clusters
luad_hier_stat_2 <- rbind(
  append('luad_com_hier_2', summary(as.factor(luad_com2_hier))),
  append('luad_com_hier_knn_2', summary(as.factor(luad_cluster2_com_knn_hier))),
  append('luad_com_hier_lls_2', summary(as.factor(luad_cluster2_com_lls_hier))),
  append('luad_com_hier_lsa_2', summary(as.factor(luad_cluster2_com_lsa_hier))),
  append('luad_com_hier_rand_2', summary(as.factor(luad_cluster2_com_rand_hier))),
  append('luad_union_hier_knn_2', summary(as.factor(luad_cluster2_union_knn_hier))),
  append('luad_union_hier_lls_2', summary(as.factor(luad_cluster2_union_lls_hier))),
  append('luad_union_hier_lsa_2', summary(as.factor(luad_cluster2_union_lsa_hier))),
  append('luad_union_hier_rand_2', summary(as.factor(luad_cluster2_union_rand_hier)))
)


# icluster############################################################################################

# icluster 5 clusters
luad_iclust_stat_5 <- rbind(
  append('luad_com_iclust_5', summary(as.factor(luad_com5_iclust))),
  append('luad_com_iclust_knn_5', summary(as.factor(luad_cluster5_com_knn_iclust))),
  append('luad_com_iclust_lls_5', summary(as.factor(luad_cluster5_com_lls_iclust))),
  append('luad_com_iclust_lsa_5', summary(as.factor(luad_cluster5_com_lsa_iclust))),
  append('luad_com_iclust_rand_5', summary(as.factor(luad_cluster5_com_rand_iclust))),
  append('luad_union_iclust_knn_5', summary(as.factor(luad_cluster5_union_knn_iclust))),
  append('luad_union_iclust_lls_5', summary(as.factor(luad_cluster5_union_lls_iclust))),
  append('luad_union_iclust_lsa_5', summary(as.factor(luad_cluster5_union_lsa_iclust))),
  append('luad_union_iclust_rand_5', summary(as.factor(luad_cluster5_union_rand_iclust)))
)

# icluster 4 clusters
luad_iclust_stat_4 <- rbind(
  append('luad_com_iclust_4', summary(as.factor(luad_com4_iclust))),
  append('luad_com_iclust_knn_4', summary(as.factor(luad_cluster4_com_knn_iclust))),
  append('luad_com_iclust_lls_4', summary(as.factor(luad_cluster4_com_lls_iclust))),
  append('luad_com_iclust_lsa_4', summary(as.factor(luad_cluster4_com_lsa_iclust))),
  append('luad_com_iclust_rand_4', summary(as.factor(luad_cluster4_com_rand_iclust))),
  append('luad_union_iclust_knn_4', summary(as.factor(luad_cluster4_union_knn_iclust))),
  append('luad_union_iclust_lls_4', summary(as.factor(luad_cluster4_union_lls_iclust))),
  append('luad_union_iclust_lsa_4', summary(as.factor(luad_cluster4_union_lsa_iclust))),
  append('luad_union_iclust_rand_4', summary(as.factor(luad_cluster4_union_rand_iclust)))
)

# icluster 3 clusters
luad_iclust_stat_3 <- rbind(
  append('luad_com_iclust_3', summary(as.factor(luad_com3_iclust))),
  append('luad_com_iclust_knn_3', summary(as.factor(luad_cluster3_com_knn_iclust))),
  append('luad_com_iclust_lls_3', summary(as.factor(luad_cluster3_com_lls_iclust))),
  append('luad_com_iclust_lsa_3', summary(as.factor(luad_cluster3_com_lsa_iclust))),
  append('luad_com_iclust_rand_3', summary(as.factor(luad_cluster3_com_rand_iclust))),
  append('luad_union_iclust_knn_3', summary(as.factor(luad_cluster3_union_knn_iclust))),
  append('luad_union_iclust_lls_3', summary(as.factor(luad_cluster3_union_lls_iclust))),
  append('luad_union_iclust_lsa_3', summary(as.factor(luad_cluster3_union_lsa_iclust))),
  append('luad_union_iclust_rand_3', summary(as.factor(luad_cluster3_union_rand_iclust)))
)

# icluster 2 clusters
luad_iclust_stat_2 <- rbind(
  append('luad_com_iclust_2', summary(as.factor(luad_com2_iclust))),
  append('luad_com_iclust_knn_2', summary(as.factor(luad_cluster2_com_knn_iclust))),
  append('luad_com_iclust_lls_2', summary(as.factor(luad_cluster2_com_lls_iclust))),
  append('luad_com_iclust_lsa_2', summary(as.factor(luad_cluster2_com_lsa_iclust))),
  append('luad_com_iclust_rand_2', summary(as.factor(luad_cluster2_com_rand_iclust))),
  append('luad_union_iclust_knn_2', summary(as.factor(luad_cluster2_union_knn_iclust))),
  append('luad_union_iclust_lls_2', summary(as.factor(luad_cluster2_union_lls_iclust))),
  append('luad_union_iclust_lsa_2', summary(as.factor(luad_cluster2_union_lsa_iclust))),
  append('luad_union_iclust_rand_2', summary(as.factor(luad_cluster2_union_rand_iclust)))
)


# snf############################################################################################

# snf 5 clusters
luad_snf_stat_5 <- rbind(
  
  append('luad_com_snf_5', summary(as.factor(luad_com5_snf))),
  append('luad_com_snf_knn_5', summary(as.factor(luad_cluster5_com_knn_snf))),
  append('luad_com_snf_lls_5', summary(as.factor(luad_cluster5_com_lls_snf))),
  append('luad_com_snf_lsa_5', summary(as.factor(luad_cluster5_com_lsa_snf))),
  append('luad_com_snf_rand_5', summary(as.factor(luad_cluster5_com_rand_snf))),
  append('luad_com_snf_reg_5', summary(as.factor(luad_cluster5_com_reg_snf))),
  append('luad_com_snf_med_5', summary(as.factor(luad_cluster5_com_med_snf))),
  append('luad_com_snf_self_5', summary(as.factor(luad_cluster5_com_self_snf))),
  append('luad_union_snf_knn_5', summary(as.factor(luad_cluster5_union_knn_snf))),
  append('luad_union_snf_lls_5', summary(as.factor(luad_cluster5_union_lls_snf))),
  append('luad_union_snf_lsa_5', summary(as.factor(luad_cluster5_union_lsa_snf))),
  append('luad_union_snf_rand_5', summary(as.factor(luad_cluster5_union_rand_snf))),
  append('luad_union_snf_reg_5', summary(as.factor(luad_cluster5_union_reg_snf))),
  append('luad_union_snf_med_5', summary(as.factor(luad_cluster5_union_med_snf))),
  append('luad_union_snf_self_5', summary(as.factor(luad_cluster5_union_self_snf)))
  
)

# snf 4 clusters
luad_snf_stat_4 <- rbind(
  
  append('luad_com_snf_4', summary(as.factor(luad_com4_snf))),
  append('luad_com_snf_knn_4', summary(as.factor(luad_cluster4_com_knn_snf))),
  append('luad_com_snf_lls_4', summary(as.factor(luad_cluster4_com_lls_snf))),
  append('luad_com_snf_lsa_4', summary(as.factor(luad_cluster4_com_lsa_snf))),
  append('luad_com_snf_rand_4', summary(as.factor(luad_cluster4_com_rand_snf))),
  append('luad_com_snf_reg_4', summary(as.factor(luad_cluster4_com_reg_snf))),
  append('luad_com_snf_med_4', summary(as.factor(luad_cluster4_com_med_snf))),
  append('luad_com_snf_self_4', summary(as.factor(luad_cluster4_com_self_snf))),
  append('luad_union_snf_knn_4', summary(as.factor(luad_cluster4_union_knn_snf))),
  append('luad_union_snf_lls_4', summary(as.factor(luad_cluster4_union_lls_snf))),
  append('luad_union_snf_lsa_4', summary(as.factor(luad_cluster4_union_lsa_snf))),
  append('luad_union_snf_rand_4', summary(as.factor(luad_cluster4_union_rand_snf))),
  append('luad_union_snf_reg_4', summary(as.factor(luad_cluster4_union_reg_snf))),
  append('luad_union_snf_med_4', summary(as.factor(luad_cluster4_union_med_snf))),
  append('luad_union_snf_self_4', summary(as.factor(luad_cluster4_union_self_snf)))
  
)


# snf 3 clusters
luad_snf_stat_3 <- rbind(
  
  append('luad_com_snf_3', summary(as.factor(luad_com3_snf))),
  append('luad_com_snf_knn_3', summary(as.factor(luad_cluster3_com_knn_snf))),
  append('luad_com_snf_lls_3', summary(as.factor(luad_cluster3_com_lls_snf))),
  append('luad_com_snf_lsa_3', summary(as.factor(luad_cluster3_com_lsa_snf))),
  append('luad_com_snf_rand_3', summary(as.factor(luad_cluster3_com_rand_snf))),
  append('luad_com_snf_reg_3', summary(as.factor(luad_cluster3_com_reg_snf))),
  append('luad_com_snf_med_3', summary(as.factor(luad_cluster3_com_med_snf))),
  append('luad_com_snf_self_3', summary(as.factor(luad_cluster3_com_self_snf))),
  append('luad_union_snf_knn_3', summary(as.factor(luad_cluster3_union_knn_snf))),
  append('luad_union_snf_lls_3', summary(as.factor(luad_cluster3_union_lls_snf))),
  append('luad_union_snf_lsa_3', summary(as.factor(luad_cluster3_union_lsa_snf))),
  append('luad_union_snf_rand_3', summary(as.factor(luad_cluster3_union_rand_snf))),
  append('luad_union_snf_reg_3', summary(as.factor(luad_cluster3_union_reg_snf))),
  append('luad_union_snf_med_3', summary(as.factor(luad_cluster3_union_med_snf))),
  append('luad_union_snf_self_3', summary(as.factor(luad_cluster3_union_self_snf)))
  
)

# snf 2 clusters
luad_snf_stat_2 <- rbind(
  
  append('luad_com_snf_2', summary(as.factor(luad_com2_snf))),
  append('luad_com_snf_knn_2', summary(as.factor(luad_cluster2_com_knn_snf))),
  append('luad_com_snf_lls_2', summary(as.factor(luad_cluster2_com_lls_snf))),
  append('luad_com_snf_lsa_2', summary(as.factor(luad_cluster2_com_lsa_snf))),
  append('luad_com_snf_rand_2', summary(as.factor(luad_cluster2_com_rand_snf))),
  append('luad_com_snf_reg_2', summary(as.factor(luad_cluster2_com_reg_snf))),
  append('luad_com_snf_med_2', summary(as.factor(luad_cluster2_com_med_snf))),
  append('luad_com_snf_self_2', summary(as.factor(luad_cluster2_com_self_snf))),
  append('luad_union_snf_knn_2', summary(as.factor(luad_cluster2_union_knn_snf))),
  append('luad_union_snf_lls_2', summary(as.factor(luad_cluster2_union_lls_snf))),
  append('luad_union_snf_lsa_2', summary(as.factor(luad_cluster2_union_lsa_snf))),
  append('luad_union_snf_rand_2', summary(as.factor(luad_cluster2_union_rand_snf))),
  append('luad_union_snf_reg_2', summary(as.factor(luad_cluster2_union_reg_snf))),
  append('luad_union_snf_med_2', summary(as.factor(luad_cluster2_union_med_snf))),
  append('luad_union_snf_self_2', summary(as.factor(luad_cluster2_union_self_snf)))
  
)

#################################################################################################

##########################################################################################################
# lusc, compare clusters from complete, compete imputed, and union

# get lusc IDs
lusc_complete_ids <- lusc_ids[[1]]
lusc_union_ids <- lusc_ids[[2]]

# Get labels for each clustering type, from the intersection for lusc
lusc_com5_hier <- as.factor(t(read.table(paste0(completeFolder, '/5_1_1.txt'))))
lusc_com5_iclust <- as.factor(t(read.table(paste0(completeFolder, '/5_1_2.txt'))))
lusc_com5_snf <- as.factor(t(read.table(paste0(completeFolder, '/5_1_3.txt'))))

lusc_com4_hier <- as.factor(t(read.table(paste0(completeFolder, '/5_2_1.txt'))))
lusc_com4_iclust <- as.factor(t(read.table(paste0(completeFolder, '/5_2_2.txt'))))
lusc_com4_snf <- as.factor(t(read.table(paste0(completeFolder, '/5_2_3.txt'))))

lusc_com3_hier <- as.factor(t(read.table(paste0(completeFolder, '/5_3_1.txt'))))
lusc_com3_iclust <- as.factor(t(read.table(paste0(completeFolder, '/5_3_2.txt'))))
lusc_com3_snf <- as.factor(t(read.table(paste0(completeFolder, '/5_3_3.txt'))))

lusc_com2_hier <- as.factor(t(read.table(paste0(completeFolder, '/5_4_1.txt'))))
lusc_com2_iclust <- as.factor(t(read.table(paste0(completeFolder, '/5_4_2.txt'))))
lusc_com2_snf <- as.factor(t(read.table(paste0(completeFolder, '/5_4_3.txt'))))
############################################################################################3
# get labels for lusc intersection imputed
# knn cluster 1 - 4
lusc_cluster5_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster1_knn_hier.txt'))))
lusc_cluster4_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster2_knn_hier.txt'))))
lusc_cluster3_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster3_knn_hier.txt'))))
lusc_cluster2_com_knn_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster4_knn_hier.txt'))))

lusc_cluster5_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster1_knn_iclust.txt'))))
lusc_cluster4_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster2_knn_iclust.txt'))))
lusc_cluster3_com_knn_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster3_knn_iclust.txt'))))
lusc_cluster2_com_knn_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster4_knn_iclust.txt'))))

lusc_cluster5_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster1_knn_snf.txt'))))
lusc_cluster4_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster2_knn_snf.txt'))))
lusc_cluster3_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster3_knn_snf.txt'))))
lusc_cluster2_com_knn_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
lusc_cluster5_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster1_lls_hier.txt'))))
lusc_cluster4_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster2_lls_hier.txt'))))
lusc_cluster3_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster3_lls_hier.txt'))))
lusc_cluster2_com_lls_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster4_lls_hier.txt'))))

lusc_cluster5_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster1_lls_iclust.txt'))))
lusc_cluster4_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster2_lls_iclust.txt'))))
lusc_cluster3_com_lls_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster3_lls_iclust.txt'))))
lusc_cluster2_com_lls_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster4_lls_iclust.txt'))))

lusc_cluster5_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster1_lls_snf.txt'))))
lusc_cluster4_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster2_lls_snf.txt'))))
lusc_cluster3_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster3_lls_snf.txt'))))
lusc_cluster2_com_lls_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
lusc_cluster5_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster1_lsa_hier.txt'))))
lusc_cluster4_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster2_lsa_hier.txt'))))
lusc_cluster3_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster3_lsa_hier.txt'))))
lusc_cluster2_com_lsa_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster4_lsa_hier.txt'))))

lusc_cluster5_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster1_lsa_iclust.txt'))))
lusc_cluster4_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster2_lsa_iclust.txt'))))
lusc_cluster3_com_lsa_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster3_lsa_iclust.txt'))))
lusc_cluster2_com_lsa_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster4_lsa_iclust.txt'))))

lusc_cluster5_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster1_lsa_snf.txt'))))
lusc_cluster4_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster2_lsa_snf.txt'))))
lusc_cluster3_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster3_lsa_snf.txt'))))
lusc_cluster2_com_lsa_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
lusc_cluster5_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster1_rand_hier.txt'))))
lusc_cluster4_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster2_rand_hier.txt'))))
lusc_cluster3_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster3_rand_hier.txt'))))
lusc_cluster2_com_rand_hier <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster4_rand_hier.txt'))))

lusc_cluster5_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster1_rand_iclust.txt'))))
lusc_cluster4_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster2_rand_iclust.txt'))))
lusc_cluster3_com_rand_iclust<- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster3_rand_iclust.txt'))))
lusc_cluster2_com_rand_iclust <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster4_rand_iclust.txt'))))

lusc_cluster5_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster1_rand_snf.txt'))))
lusc_cluster4_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster2_rand_snf.txt'))))
lusc_cluster3_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster3_rand_snf.txt'))))
lusc_cluster2_com_rand_snf <- as.factor(t(read.table(paste0(imputeFolder, '/lusc/rand', '/cluster4_rand_snf.txt'))))


# similarity
lusc_cluster5_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster1_med_snf.txt'))))
lusc_cluster4_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster2_med_snf.txt'))))
lusc_cluster3_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster3_med_snf.txt'))))
lusc_cluster2_com_med_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster4_med_snf.txt'))))


lusc_cluster5_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster1_regress_snf.txt'))))
lusc_cluster4_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster2_regress_snf.txt'))))
lusc_cluster3_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster3_regress_snf.txt'))))
lusc_cluster2_com_reg_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster4_regress_snf.txt'))))

lusc_cluster5_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster1_self_snf.txt'))))
lusc_cluster4_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster2_self_snf.txt'))))
lusc_cluster3_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster3_self_snf.txt'))))
lusc_cluster2_com_self_snf <- as.factor(t(read.table(paste0(similarityFolder, '/lusc', '/cluster4_self_snf.txt'))))

#########################################################################################################
# Get labels for union

# knn cluster 1 - 4
lusc_cluster5_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster1_knn_hier.txt'))))
lusc_cluster4_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster2_knn_hier.txt'))))
lusc_cluster3_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster3_knn_hier.txt'))))
lusc_cluster2_union_knn_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster4_knn_hier.txt'))))

lusc_cluster5_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster1_knn_iclust.txt'))))
lusc_cluster4_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster2_knn_iclust.txt'))))
lusc_cluster3_union_knn_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster3_knn_iclust.txt'))))
lusc_cluster2_union_knn_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster4_knn_iclust.txt'))))

lusc_cluster5_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster1_knn_snf.txt'))))
lusc_cluster4_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster2_knn_snf.txt'))))
lusc_cluster3_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster3_knn_snf.txt'))))
lusc_cluster2_union_knn_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/knn', '/cluster4_knn_snf.txt'))))

# lls cluster 1 - 4
lusc_cluster5_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster1_lls_hier.txt'))))
lusc_cluster4_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster2_lls_hier.txt'))))
lusc_cluster3_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster3_lls_hier.txt'))))
lusc_cluster2_union_lls_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster4_lls_hier.txt'))))

lusc_cluster5_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster1_lls_iclust.txt'))))
lusc_cluster4_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster2_lls_iclust.txt'))))
lusc_cluster3_union_lls_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster3_lls_iclust.txt'))))
lusc_cluster2_union_lls_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster4_lls_iclust.txt'))))

lusc_cluster5_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster1_lls_snf.txt'))))
lusc_cluster4_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster2_lls_snf.txt'))))
lusc_cluster3_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster3_lls_snf.txt'))))
lusc_cluster2_union_lls_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lls', '/cluster4_lls_snf.txt'))))

# lsa cluster 1 - 4
lusc_cluster5_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster1_lsa_hier.txt'))))
lusc_cluster4_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster2_lsa_hier.txt'))))
lusc_cluster3_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster3_lsa_hier.txt'))))
lusc_cluster2_union_lsa_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster4_lsa_hier.txt'))))

lusc_cluster5_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster1_lsa_iclust.txt'))))
lusc_cluster4_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster2_lsa_iclust.txt'))))
lusc_cluster3_union_lsa_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster3_lsa_iclust.txt'))))
lusc_cluster2_union_lsa_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster4_lsa_iclust.txt'))))

lusc_cluster5_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster1_lsa_snf.txt'))))
lusc_cluster4_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster2_lsa_snf.txt'))))
lusc_cluster3_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster3_lsa_snf.txt'))))
lusc_cluster2_union_lsa_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/lsa', '/cluster4_lsa_snf.txt'))))

# rand cluster 1 - 4
lusc_cluster5_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster1_rand_hier.txt'))))
lusc_cluster4_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster2_rand_hier.txt'))))
lusc_cluster3_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster3_rand_hier.txt'))))
lusc_cluster2_union_rand_hier <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster4_rand_hier.txt'))))

lusc_cluster5_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster1_rand_iclust.txt'))))
lusc_cluster4_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster2_rand_iclust.txt'))))
lusc_cluster3_union_rand_iclust<- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster3_rand_iclust.txt'))))
lusc_cluster2_union_rand_iclust <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster4_rand_iclust.txt'))))

lusc_cluster5_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster1_rand_snf.txt'))))
lusc_cluster4_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster2_rand_snf.txt'))))
lusc_cluster3_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster3_rand_snf.txt'))))
lusc_cluster2_union_rand_snf <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lusc/rand', '/cluster4_rand_snf.txt'))))


# similarity
lusc_cluster5_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster1_med_snf.txt'))))
lusc_cluster4_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster2_med_snf.txt'))))
lusc_cluster3_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster3_med_snf.txt'))))
lusc_cluster2_union_med_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster4_med_snf.txt'))))


lusc_cluster5_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster1_regress_snf.txt'))))
lusc_cluster4_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster2_regress_snf.txt'))))
lusc_cluster3_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster3_regress_snf.txt'))))
lusc_cluster2_union_reg_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster4_regress_snf.txt'))))

lusc_cluster5_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster1_self_snf.txt'))))
lusc_cluster4_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster2_self_snf.txt'))))
lusc_cluster3_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster3_self_snf.txt'))))
lusc_cluster2_union_self_snf <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lusc', '/cluster4_self_snf.txt'))))


##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical############################################################################################

# Hierarchical 5 clusters
lusc_hier_stat_5 <- rbind(
  append('lusc_com_hier_5', summary(as.factor(lusc_com5_hier))),
  append('lusc_com_hier_knn_5', summary(as.factor(lusc_cluster5_com_knn_hier))),
  append('lusc_com_hier_lls_5', summary(as.factor(lusc_cluster5_com_lls_hier))),
  append('lusc_com_hier_lsa_5', summary(as.factor(lusc_cluster5_com_lsa_hier))),
  append('lusc_com_hier_rand_5', summary(as.factor(lusc_cluster5_com_rand_hier))),
  append('lusc_union_hier_knn_5', summary(as.factor(lusc_cluster5_union_knn_hier))),
  append('lusc_union_hier_lls_5', summary(as.factor(lusc_cluster5_union_lls_hier))),
  append('lusc_union_hier_lsa_5', summary(as.factor(lusc_cluster5_union_lsa_hier))),
  append('lusc_union_hier_rand_5', summary(as.factor(lusc_cluster5_union_rand_hier)))
)

# Hierarchical 4 clusters
lusc_hier_stat_4 <- rbind(
  append('lusc_com_hier_4', summary(as.factor(lusc_com4_hier))),
  append('lusc_com_hier_knn_4', summary(as.factor(lusc_cluster4_com_knn_hier))),
  append('lusc_com_hier_lls_4', summary(as.factor(lusc_cluster4_com_lls_hier))),
  append('lusc_com_hier_lsa_4', summary(as.factor(lusc_cluster4_com_lsa_hier))),
  append('lusc_com_hier_rand_4', summary(as.factor(lusc_cluster4_com_rand_hier))),
  append('lusc_union_hier_knn_4', summary(as.factor(lusc_cluster4_union_knn_hier))),
  append('lusc_union_hier_lls_4', summary(as.factor(lusc_cluster4_union_lls_hier))),
  append('lusc_union_hier_lsa_4', summary(as.factor(lusc_cluster4_union_lsa_hier))),
  append('lusc_union_hier_rand_4', summary(as.factor(lusc_cluster4_union_rand_hier)))
)

# Hierarchical 3 clusters
lusc_hier_stat_3 <- rbind(
  append('lusc_com_hier_3', summary(as.factor(lusc_com3_hier))),
  append('lusc_com_hier_knn_3', summary(as.factor(lusc_cluster3_com_knn_hier))),
  append('lusc_com_hier_lls_3', summary(as.factor(lusc_cluster3_com_lls_hier))),
  append('lusc_com_hier_lsa_3', summary(as.factor(lusc_cluster3_com_lsa_hier))),
  append('lusc_com_hier_rand_3', summary(as.factor(lusc_cluster3_com_rand_hier))),
  append('lusc_union_hier_knn_3', summary(as.factor(lusc_cluster3_union_knn_hier))),
  append('lusc_union_hier_lls_3', summary(as.factor(lusc_cluster3_union_lls_hier))),
  append('lusc_union_hier_lsa_3', summary(as.factor(lusc_cluster3_union_lsa_hier))),
  append('lusc_union_hier_rand_3', summary(as.factor(lusc_cluster3_union_rand_hier)))
)

# Hierarchical 2 clusters
lusc_hier_stat_2 <- rbind(
  append('lusc_com_hier_2', summary(as.factor(lusc_com2_hier))),
  append('lusc_com_hier_knn_2', summary(as.factor(lusc_cluster2_com_knn_hier))),
  append('lusc_com_hier_lls_2', summary(as.factor(lusc_cluster2_com_lls_hier))),
  append('lusc_com_hier_lsa_2', summary(as.factor(lusc_cluster2_com_lsa_hier))),
  append('lusc_com_hier_rand_2', summary(as.factor(lusc_cluster2_com_rand_hier))),
  append('lusc_union_hier_knn_2', summary(as.factor(lusc_cluster2_union_knn_hier))),
  append('lusc_union_hier_lls_2', summary(as.factor(lusc_cluster2_union_lls_hier))),
  append('lusc_union_hier_lsa_2', summary(as.factor(lusc_cluster2_union_lsa_hier))),
  append('lusc_union_hier_rand_2', summary(as.factor(lusc_cluster2_union_rand_hier)))
)


# icluster############################################################################################

# icluster 5 clusters
lusc_iclust_stat_5 <- rbind(
  append('lusc_com_iclust_5', summary(as.factor(lusc_com5_iclust))),
  append('lusc_com_iclust_knn_5', summary(as.factor(lusc_cluster5_com_knn_iclust))),
  append('lusc_com_iclust_lls_5', summary(as.factor(lusc_cluster5_com_lls_iclust))),
  append('lusc_com_iclust_lsa_5', summary(as.factor(lusc_cluster5_com_lsa_iclust))),
  append('lusc_com_iclust_rand_5', summary(as.factor(lusc_cluster5_com_rand_iclust))),
  append('lusc_union_iclust_knn_5', summary(as.factor(lusc_cluster5_union_knn_iclust))),
  append('lusc_union_iclust_lls_5', summary(as.factor(lusc_cluster5_union_lls_iclust))),
  append('lusc_union_iclust_lsa_5', summary(as.factor(lusc_cluster5_union_lsa_iclust))),
  append('lusc_union_iclust_rand_5', summary(as.factor(lusc_cluster5_union_rand_iclust)))
)

# icluster 4 clusters
lusc_iclust_stat_4 <- rbind(
  append('lusc_com_iclust_4', summary(as.factor(lusc_com4_iclust))),
  append('lusc_com_iclust_knn_4', summary(as.factor(lusc_cluster4_com_knn_iclust))),
  append('lusc_com_iclust_lls_4', summary(as.factor(lusc_cluster4_com_lls_iclust))),
  append('lusc_com_iclust_lsa_4', summary(as.factor(lusc_cluster4_com_lsa_iclust))),
  append('lusc_com_iclust_rand_4', summary(as.factor(lusc_cluster4_com_rand_iclust))),
  append('lusc_union_iclust_knn_4', summary(as.factor(lusc_cluster4_union_knn_iclust))),
  append('lusc_union_iclust_lls_4', summary(as.factor(lusc_cluster4_union_lls_iclust))),
  append('lusc_union_iclust_lsa_4', summary(as.factor(lusc_cluster4_union_lsa_iclust))),
  append('lusc_union_iclust_rand_4', summary(as.factor(lusc_cluster4_union_rand_iclust)))
)

# icluster 3 clusters
lusc_iclust_stat_3 <- rbind(
  append('lusc_com_iclust_3', summary(as.factor(lusc_com3_iclust))),
  append('lusc_com_iclust_knn_3', summary(as.factor(lusc_cluster3_com_knn_iclust))),
  append('lusc_com_iclust_lls_3', summary(as.factor(lusc_cluster3_com_lls_iclust))),
  append('lusc_com_iclust_lsa_3', summary(as.factor(lusc_cluster3_com_lsa_iclust))),
  append('lusc_com_iclust_rand_3', summary(as.factor(lusc_cluster3_com_rand_iclust))),
  append('lusc_union_iclust_knn_3', summary(as.factor(lusc_cluster3_union_knn_iclust))),
  append('lusc_union_iclust_lls_3', summary(as.factor(lusc_cluster3_union_lls_iclust))),
  append('lusc_union_iclust_lsa_3', summary(as.factor(lusc_cluster3_union_lsa_iclust))),
  append('lusc_union_iclust_rand_3', summary(as.factor(lusc_cluster3_union_rand_iclust)))
)

# icluster 2 clusters
lusc_iclust_stat_2 <- rbind(
  append('lusc_com_iclust_2', summary(as.factor(lusc_com2_iclust))),
  append('lusc_com_iclust_knn_2', summary(as.factor(lusc_cluster2_com_knn_iclust))),
  append('lusc_com_iclust_lls_2', summary(as.factor(lusc_cluster2_com_lls_iclust))),
  append('lusc_com_iclust_lsa_2', summary(as.factor(lusc_cluster2_com_lsa_iclust))),
  append('lusc_com_iclust_rand_2', summary(as.factor(lusc_cluster2_com_rand_iclust))),
  append('lusc_union_iclust_knn_2', summary(as.factor(lusc_cluster2_union_knn_iclust))),
  append('lusc_union_iclust_lls_2', summary(as.factor(lusc_cluster2_union_lls_iclust))),
  append('lusc_union_iclust_lsa_2', summary(as.factor(lusc_cluster2_union_lsa_iclust))),
  append('lusc_union_iclust_rand_2', summary(as.factor(lusc_cluster2_union_rand_iclust)))
)


# snf############################################################################################

# snf 5 clusters
lusc_snf_stat_5 <- rbind(
  
  append('lusc_com_snf_5', summary(as.factor(lusc_com5_snf))),
  append('lusc_com_snf_knn_5', summary(as.factor(lusc_cluster5_com_knn_snf))),
  append('lusc_com_snf_lls_5', summary(as.factor(lusc_cluster5_com_lls_snf))),
  append('lusc_com_snf_lsa_5', summary(as.factor(lusc_cluster5_com_lsa_snf))),
  append('lusc_com_snf_rand_5', summary(as.factor(lusc_cluster5_com_rand_snf))),
  append('lusc_com_snf_reg_5', summary(as.factor(lusc_cluster5_com_reg_snf))),
  append('lusc_com_snf_med_5', summary(as.factor(lusc_cluster5_com_med_snf))),
  append('lusc_com_snf_self_5', summary(as.factor(lusc_cluster5_com_self_snf))),
  append('lusc_union_snf_knn_5', summary(as.factor(lusc_cluster5_union_knn_snf))),
  append('lusc_union_snf_lls_5', summary(as.factor(lusc_cluster5_union_lls_snf))),
  append('lusc_union_snf_lsa_5', summary(as.factor(lusc_cluster5_union_lsa_snf))),
  append('lusc_union_snf_rand_5', summary(as.factor(lusc_cluster5_union_rand_snf))),
  append('lusc_union_snf_reg_5', summary(as.factor(lusc_cluster5_union_reg_snf))),
  append('lusc_union_snf_med_5', summary(as.factor(lusc_cluster5_union_med_snf))),
  append('lusc_union_snf_self_5', summary(as.factor(lusc_cluster5_union_self_snf)))
  
)

# snf 4 clusters
lusc_snf_stat_4 <- rbind(
  
  append('lusc_com_snf_4', summary(as.factor(lusc_com4_snf))),
  append('lusc_com_snf_knn_4', summary(as.factor(lusc_cluster4_com_knn_snf))),
  append('lusc_com_snf_lls_4', summary(as.factor(lusc_cluster4_com_lls_snf))),
  append('lusc_com_snf_lsa_4', summary(as.factor(lusc_cluster4_com_lsa_snf))),
  append('lusc_com_snf_rand_4', summary(as.factor(lusc_cluster4_com_rand_snf))),
  append('lusc_com_snf_reg_4', summary(as.factor(lusc_cluster4_com_reg_snf))),
  append('lusc_com_snf_med_4', summary(as.factor(lusc_cluster4_com_med_snf))),
  append('lusc_com_snf_self_4', summary(as.factor(lusc_cluster4_com_self_snf))),
  append('lusc_union_snf_knn_4', summary(as.factor(lusc_cluster4_union_knn_snf))),
  append('lusc_union_snf_lls_4', summary(as.factor(lusc_cluster4_union_lls_snf))),
  append('lusc_union_snf_lsa_4', summary(as.factor(lusc_cluster4_union_lsa_snf))),
  append('lusc_union_snf_rand_4', summary(as.factor(lusc_cluster4_union_rand_snf))),
  append('lusc_union_snf_reg_4', summary(as.factor(lusc_cluster4_union_reg_snf))),
  append('lusc_union_snf_med_4', summary(as.factor(lusc_cluster4_union_med_snf))),
  append('lusc_union_snf_self_4', summary(as.factor(lusc_cluster4_union_self_snf)))
  
)


# snf 3 clusters
lusc_snf_stat_3 <- rbind(
  
  append('lusc_com_snf_3', summary(as.factor(lusc_com3_snf))),
  append('lusc_com_snf_knn_3', summary(as.factor(lusc_cluster3_com_knn_snf))),
  append('lusc_com_snf_lls_3', summary(as.factor(lusc_cluster3_com_lls_snf))),
  append('lusc_com_snf_lsa_3', summary(as.factor(lusc_cluster3_com_lsa_snf))),
  append('lusc_com_snf_rand_3', summary(as.factor(lusc_cluster3_com_rand_snf))),
  append('lusc_com_snf_reg_3', summary(as.factor(lusc_cluster3_com_reg_snf))),
  append('lusc_com_snf_med_3', summary(as.factor(lusc_cluster3_com_med_snf))),
  append('lusc_com_snf_self_3', summary(as.factor(lusc_cluster3_com_self_snf))),
  append('lusc_union_snf_knn_3', summary(as.factor(lusc_cluster3_union_knn_snf))),
  append('lusc_union_snf_lls_3', summary(as.factor(lusc_cluster3_union_lls_snf))),
  append('lusc_union_snf_lsa_3', summary(as.factor(lusc_cluster3_union_lsa_snf))),
  append('lusc_union_snf_rand_3', summary(as.factor(lusc_cluster3_union_rand_snf))),
  append('lusc_union_snf_reg_3', summary(as.factor(lusc_cluster3_union_reg_snf))),
  append('lusc_union_snf_med_3', summary(as.factor(lusc_cluster3_union_med_snf))),
  append('lusc_union_snf_self_3', summary(as.factor(lusc_cluster3_union_self_snf)))
  
)

# snf 2 clusters
lusc_snf_stat_2 <- rbind(
  
  append('lusc_com_snf_2', summary(as.factor(lusc_com2_snf))),
  append('lusc_com_snf_knn_2', summary(as.factor(lusc_cluster2_com_knn_snf))),
  append('lusc_com_snf_lls_2', summary(as.factor(lusc_cluster2_com_lls_snf))),
  append('lusc_com_snf_lsa_2', summary(as.factor(lusc_cluster2_com_lsa_snf))),
  append('lusc_com_snf_rand_2', summary(as.factor(lusc_cluster2_com_rand_snf))),
  append('lusc_com_snf_reg_2', summary(as.factor(lusc_cluster2_com_reg_snf))),
  append('lusc_com_snf_med_2', summary(as.factor(lusc_cluster2_com_med_snf))),
  append('lusc_com_snf_self_2', summary(as.factor(lusc_cluster2_com_self_snf))),
  append('lusc_union_snf_knn_2', summary(as.factor(lusc_cluster2_union_knn_snf))),
  append('lusc_union_snf_lls_2', summary(as.factor(lusc_cluster2_union_lls_snf))),
  append('lusc_union_snf_lsa_2', summary(as.factor(lusc_cluster2_union_lsa_snf))),
  append('lusc_union_snf_rand_2', summary(as.factor(lusc_cluster2_union_rand_snf))),
  append('lusc_union_snf_reg_2', summary(as.factor(lusc_cluster2_union_reg_snf))),
  append('lusc_union_snf_med_2', summary(as.factor(lusc_cluster2_union_med_snf))),
  append('lusc_union_snf_self_2', summary(as.factor(lusc_cluster2_union_self_snf)))
  
)

#################################################################################################
