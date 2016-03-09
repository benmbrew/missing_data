#########################################################################################################
# This script will compare labels from the intersection with those from the union. 

# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)
library(dplyr)

#########################################################################################################
# initiate folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
completeFolder <- paste(projectFolder, "Scripts",
                        "06_Two_Thousand_Features",
                        "cluster_complete_data", 'Results', 'Labels',sep="/")
imputeFolder <- paste(projectFolder, "Scripts",
                    "06_Int_Labels",
                    "evaluate_original_imputation", sep="/")
similarityFolder <- paste(projectFolder, "Scripts",
                      "06_Int_Labels",
                      "evaluate_original_similarity", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Cluster_Compare", sep="/")
idsFolder <- paste(testFolder, "ids", sep="/")

#########################################################################################################
# Get ids from data for each cancer
source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))

# Load the original data
loadIDs <- function(cancer){
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
  

brca_ids <-loadIDs(cancer = 'BRCA')
kirc_ids <-loadIDs(cancer = 'KIRC')
lihc_ids <-loadIDs(cancer = 'LIHC')
luad_ids <-loadIDs(cancer = 'LUAD')

# unlist ids into complete and union
brca_complete_ids <- brca_ids[[1]]
brca_union_ids <- brca_ids[[2]]
kirc_complete_ids <- kirc_ids[[1]]
kirc_union_ids <- kirc_ids[[2]]
lihc_complete_ids <- lihc_ids[[1]]
lihc_union_ids <- lihc_ids[[2]]
luad_complete_ids <- luad_ids[[1]]
luad_union_ids <- luad_ids[[2]]


##########################################################################################################
# Get labels for each clustering type, from the intersection for brca
brca_com_hier <- read.table(paste0(completeFolder, '/1_1.txt'))
brca_com_iclus <- read.table(paste0(completeFolder, '/1_2.txt'))
brca_com_snf <- read.table(paste0(completeFolder, '/1_3.txt'))

# Get labels for each clustering type, from the intersection for kirc
kirc_com_hier <- read.table(paste0(completeFolder, '/2_1.txt'))
kirc_com_iclus <- read.table(paste0(completeFolder, '/2_2.txt'))
kirc_com_snf <- read.table(paste0(completeFolder, '/2_3.txt'))

# Get labels for each clustering type, from the intersection for lihc
lihc_com_hier <- read.table(paste0(completeFolder, '/3_1.txt'))
lihc_com_iclus <- read.table(paste0(completeFolder, '/3_2.txt'))
lihc_com_snf <- read.table(paste0(completeFolder, '/3_3.txt'))

# Get labels for each clustering type, from the intersection for luad
luad_com_hier <- read.table(paste0(completeFolder, '/4_1.txt'))
luad_com_iclus <- read.table(paste0(completeFolder, '/4_2.txt'))
luad_com_snf <- read.table(paste0(completeFolder, '/4_3.txt'))

##########################################################################################################
# Get labels for each clustering and imputation type from the union 

# get labels for brca union
brca_union_hier_knn
brca_union_hier_lls
brca_union_hier_lsa
brca_union_hier_rand

brca_union_iclus_knn
brca_union_iclus_lls
brca_union_iclus_lsa
brca_union_iclus_rand

brca_union_snf_knn
brca_union_snf_lls
brca_union_snf_lsa
brca_union_snf_rand

brca_union_snf_self
brca_union_snf_med
brca_union_snf_reg

# get labels for kirc union
kirc_union_hier_knn
kirc_union_hier_lls
kirc_union_hier_lsa
kirc_union_hier_rand

kirc_union_iclus_knn
kirc_union_iclus_lls
kirc_union_iclus_lsa
kirc_union_iclus_rand

kirc_union_snf_knn
kirc_union_snf_lls
kirc_union_snf_lsa
kirc_union_snf_rand

kirc_union_snf_self
kirc_union_snf_med
kirc_union_snf_reg

# get labels for lihc union
lihc_union_hier_knn
lihc_union_hier_lls
lihc_union_hier_lsa
lihc_union_hier_rand

lihc_union_iclus_knn
lihc_union_iclus_lls
lihc_union_iclus_lsa
lihc_union_iclus_rand

lihc_union_snf_knn
lihc_union_snf_lls
lihc_union_snf_lsa
lihc_union_snf_rand

lihc_union_snf_self
lihc_union_snf_med
lihc_union_snf_reg

# get labels for luad union
luad_union_hier_knn
luad_union_hier_lls
luad_union_hier_lsa
luad_union_hier_rand

luad_union_iclus_knn
luad_union_iclus_lls
luad_union_iclus_lsa
luad_union_iclus_rand

luad_union_snf_knn
luad_union_snf_lls
luad_union_snf_lsa
luad_union_snf_rand

luad_union_snf_self
luad_union_snf_med
luad_union_snf_reg








