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
imputeOrigFolder <- paste(projectFolder, "Scripts",
                    "06_Int_Labels",
                    "evaluate_original_imputation/Results/Clustering", sep="/")
similarityOrigFolder <- paste(projectFolder, "Scripts",
                      "06_Int_Labels",
                      "evaluate_original_similarity/Results/Similarity", sep="/")
imputeFolder <- paste(projectFolder, "Scripts",
                      "06_Int_Labels",
                      "evaluate_imputation/Results/Clustering", sep="/")
similarityFolder <- paste(projectFolder, "Scripts",
                          "06_Int_Labels",
                          "evaluate_similarity/Results/Similarity", sep="/")
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
brca_com_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_1.txt'))))
brca_com_iclus <- as.factor(t(read.table(paste0(completeFolder, '/1_2.txt'))))
brca_com_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_3.txt'))))

# get labels for brca intersection imputed
# hier
brca_com_hier_knn <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_knn_hier_comp.txt'))))
brca_com_hier_lls <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_lls_hier_comp.txt'))))
brca_com_hier_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_lsa_hier_comp.txt'))))
brca_com_hier_rand <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_rand_hier_comp.txt'))))

#iclust
brca_com_iclus_knn <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_knn_iclust_comp.txt'))))
brca_com_iclus_lls <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_lls_iclust_comp.txt'))))
brca_com_iclus_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_lsa_iclust_comp.txt'))))
brca_com_iclus_rand <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_rand_iclust_comp.txt'))))

#SNF
brca_com_snf_knn <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_knn_snf_comp.txt'))))
brca_com_snf_lls <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_lls_snf_comp.txt'))))
brca_com_snf_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_lsa_snf_comp.txt'))))
brca_com_snf_rand <- as.factor(t(read.table(paste0(imputeFolder, '/brca', '/brca_rand_snf_comp.txt'))))
brca_com_snf_self <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/brca_self_com_labels.txt'))))
brca_com_snf_med <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/brca_med_com_labels.txt' ))))
brca_com_snf_reg <- as.factor(t(read.table(paste0(similarityFolder, '/brca', '/brca_reg_com_labels.txt' ))))

###################################################################################################
# get labels for brca union

# Hier
brca_union_hier_knn  <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_knn_hier_union.txt'))))
brca_union_hier_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_lls_hier_union.txt'))))
brca_union_hier_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_lsa_hier_union.txt'))))
brca_union_hier_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_rand_hier_union.txt'))))

#iclust
brca_union_iclus_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_knn_iclust_union.txt'))))
brca_union_iclus_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_lls_iclust_union.txt'))))
brca_union_iclus_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_lsa_iclust_union.txt'))))
brca_union_iclus_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_rand_iclust_union.txt'))))


# SNF
brca_union_snf_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_knn_snf_union.txt'))))
brca_union_snf_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_lls_snf_union.txt'))))
brca_union_snf_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_lsa_snf_union.txt'))))
brca_union_snf_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/brca/', 'brca_rand_snf_union.txt'))))
brca_union_snf_self <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/brca_self_union_labels.txt' ))))
brca_union_snf_med <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/brca_med_union_labels.txt' ))))
brca_union_snf_reg <- as.factor(t(read.table(paste0(similarityOrigFolder, '/brca', '/brca_reg_union_labels.txt' ))))

##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical 
brca_hier_stat <- rbind(
  append('brca_com_hier_knn', summary(as.factor(brca_com_hier_knn))),
  append('brca_com_hier_lls', summary(as.factor(brca_com_hier_lls))),
  append('brca_com_hier_lsa', summary(as.factor(brca_com_hier_lsa))),
  append('brca_com_hier_rand', summary(as.factor(brca_com_hier_rand))),
  append('brca_union_hier_knn', summary(as.factor(brca_union_hier_knn))),
  append('brca_union_hier_lls', summary(as.factor(brca_union_hier_lls))),
  append('brca_union_hier_lsa', summary(as.factor(brca_union_hier_lsa))),
  append('brca_union_hier_rand', summary(as.factor(brca_union_hier_rand)))
  
)

# Icluster
brca_iclust_stat <- rbind(
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
  data_intersect_union <- data2[data2$id %in% data1$id,]
  data_intersect <- left_join(data_intersect_union, data1, by = 'id')
  common_clusters <- matrix(,0,5)
  cluster_compare <- matrix(,0,5)
  for(i in 1:numClus){
    sub_union <- data_intersect[data_intersect$label.x == i,]
    for(j in 1:numClus){
      sub_intersect <- data_intersect[data_intersect$label.y == j,]
      cluster_compare[j] <- round((sum(sub_union$id %in% sub_intersect$id)/nrow(sub_union))*100, 2)
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

##########################################################################################################
# kirc, compare clusters from complete, compete imputed, and union

# get kirc IDs
kirc_complete_ids <- kirc_ids[[1]]
kirc_union_ids <- kirc_ids[[2]]

# Get labels for each clustering type, from the intersection for kirc
kirc_com_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_1.txt'))))
kirc_com_iclus <- as.factor(t(read.table(paste0(completeFolder, '/1_2.txt'))))
kirc_com_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_3.txt'))))

# get labels for kirc intersection imputed
# hier
kirc_com_hier_knn <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_knn_hier_comp.txt'))))
kirc_com_hier_lls <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_lls_hier_comp.txt'))))
kirc_com_hier_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_lsa_hier_comp.txt'))))
kirc_com_hier_rand <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_rand_hier_comp.txt'))))

#iclust
kirc_com_iclus_knn <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_knn_iclust_comp.txt'))))
kirc_com_iclus_lls <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_lls_iclust_comp.txt'))))
kirc_com_iclus_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_lsa_iclust_comp.txt'))))
kirc_com_iclus_rand <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_rand_iclust_comp.txt'))))

#SNF
kirc_com_snf_knn <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_knn_snf_comp.txt'))))
kirc_com_snf_lls <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_lls_snf_comp.txt'))))
kirc_com_snf_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_lsa_snf_comp.txt'))))
kirc_com_snf_rand <- as.factor(t(read.table(paste0(imputeFolder, '/kirc', '/kirc_rand_snf_comp.txt'))))
kirc_com_snf_self <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/kirc_self_com_labels.txt'))))
kirc_com_snf_med <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/kirc_med_com_labels.txt' ))))
kirc_com_snf_reg <- as.factor(t(read.table(paste0(similarityFolder, '/kirc', '/kirc_reg_com_labels.txt' ))))

###################################################################################################
# get labels for kirc union

# Hier
kirc_union_hier_knn  <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_knn_hier_union.txt'))))
kirc_union_hier_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_lls_hier_union.txt'))))
kirc_union_hier_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_lsa_hier_union.txt'))))
kirc_union_hier_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_rand_hier_union.txt'))))

#iclust
kirc_union_iclus_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_knn_iclust_union.txt'))))
kirc_union_iclus_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_lls_iclust_union.txt'))))
kirc_union_iclus_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_lsa_iclust_union.txt'))))
kirc_union_iclus_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_rand_iclust_union.txt'))))


# SNF
kirc_union_snf_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_knn_snf_union.txt'))))
kirc_union_snf_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_lls_snf_union.txt'))))
kirc_union_snf_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_lsa_snf_union.txt'))))
kirc_union_snf_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/kirc/', 'kirc_rand_snf_union.txt'))))
kirc_union_snf_self <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/kirc_self_union_labels.txt' ))))
kirc_union_snf_med <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/kirc_med_union_labels.txt' ))))
kirc_union_snf_reg <- as.factor(t(read.table(paste0(similarityOrigFolder, '/kirc', '/kirc_reg_union_labels.txt' ))))

##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical 
kirc_hier_stat <- rbind(
  append('kirc_com_hier_knn', summary(as.factor(kirc_com_hier_knn))),
  append('kirc_com_hier_lls', summary(as.factor(kirc_com_hier_lls))),
  append('kirc_com_hier_lsa', summary(as.factor(kirc_com_hier_lsa))),
  append('kirc_com_hier_rand', summary(as.factor(kirc_com_hier_rand))),
  append('kirc_union_hier_knn', summary(as.factor(kirc_union_hier_knn))),
  append('kirc_union_hier_lls', summary(as.factor(kirc_union_hier_lls))),
  append('kirc_union_hier_lsa', summary(as.factor(kirc_union_hier_lsa))),
  append('kirc_union_hier_rand', summary(as.factor(kirc_union_hier_rand)))
  
)

# Icluster
kirc_iclust_stat <- rbind(
  append('kirc_com_iclust_knn', summary(as.factor(kirc_com_iclus_knn))),
  append('kirc_com_iclust_lls', summary(as.factor(kirc_com_iclus_lls))),
  append('kirc_com_iclust_lsa', summary(as.factor(kirc_com_iclus_lsa))),
  append('kirc_com_iclust_rand', summary(as.factor(kirc_com_iclus_rand))),
  append('kirc_union_iclust_knn', summary(as.factor(kirc_union_iclus_knn))),
  append('kirc_union_iclust_lls', summary(as.factor(kirc_union_iclus_lls))),
  append('kirc_union_iclust_lsa', summary(as.factor(kirc_union_iclus_lsa))),
  append('kirc_union_iclust_rand', summary(as.factor(kirc_union_iclus_rand)))
  
)


# SNF
kirc_snf_stat <- rbind(
  append('kirc_com_snf', summary(as.factor(kirc_com_snf))),
  append('kirc_com_snf_self', summary(as.factor(kirc_com_snf_self))),
  append('kirc_com_snf_med', summary(as.factor(kirc_com_snf_med))),
  append('kirc_com_snf_reg', summary(as.factor(kirc_com_snf_reg))),
  append('kirc_union_snf_self', summary(as.factor(kirc_union_snf_self))),
  append('kirc_union_snf_med', summary(as.factor(kirc_union_snf_med))),
  append('kirc_union_snf_reg', summary(as.factor(kirc_union_snf_reg)))
)



# Merge ids with labels and get clinical information for each cluster

# look at groups of ids in complete, complete_impute, and union
# join ids and then join each label groupd
kirc_com_hier <- as.data.frame(cbind(label = kirc_com_hier, id = kirc_complete_ids))
kirc_com_hier_knn <- as.data.frame(cbind(label = kirc_com_hier_knn, id = kirc_complete_ids))
kirc_com_hier_lls <- as.data.frame(cbind(label = kirc_com_hier_lls, id = kirc_complete_ids))
kirc_com_hier_lsa <- as.data.frame(cbind(label = kirc_com_hier_lsa, id = kirc_complete_ids))
kirc_com_hier_rand <-  as.data.frame(cbind(label = kirc_union_hier_rand, id = kirc_union_ids))
kirc_union_hier_knn <- as.data.frame(cbind(label = kirc_union_hier_knn, id = kirc_union_ids))
kirc_union_hier_lls <- as.data.frame(cbind(label = kirc_union_hier_lls, id = kirc_union_ids))
kirc_union_hier_lsa <- as.data.frame(cbind(label = kirc_union_hier_lsa, id = kirc_union_ids))
kirc_union_hier_rand <-  as.data.frame(cbind(label = kirc_union_hier_rand, id = kirc_union_ids))

kirc_com_iclus <- as.data.frame(cbind(label = kirc_com_iclus, id = kirc_complete_ids))
kirc_com_iclus_knn <- as.data.frame(cbind(label = kirc_com_iclus_knn, id = kirc_complete_ids))
kirc_com_iclus_lls <- as.data.frame(cbind(label = kirc_com_iclus_lls, id = kirc_complete_ids))
kirc_com_iclus_lsa <- as.data.frame(cbind(label = kirc_com_iclus_lsa, id = kirc_complete_ids))
kirc_com_iclus_rand <-  as.data.frame(cbind(label = kirc_union_iclus_rand, id = kirc_union_ids))
kirc_union_iclus_knn <- as.data.frame(cbind(label = kirc_union_iclus_knn, id = kirc_union_ids))
kirc_union_iclus_lls <- as.data.frame(cbind(label = kirc_union_iclus_lls, id = kirc_union_ids))
kirc_union_iclus_lsa <- as.data.frame(cbind(label = kirc_union_iclus_lsa, id = kirc_union_ids))
kirc_union_iclus_rand <-  as.data.frame(cbind(label = kirc_union_iclus_rand, id = kirc_union_ids))

kirc_com_snf <- as.data.frame(cbind(label = kirc_com_snf, id = kirc_complete_ids))
kirc_com_snf_knn <- as.data.frame(cbind(label = kirc_com_snf_knn, id = kirc_complete_ids))
kirc_com_snf_lls <- as.data.frame(cbind(label = kirc_com_snf_lls, id = kirc_complete_ids))
kirc_com_snf_lsa <- as.data.frame(cbind(label = kirc_com_snf_lsa, id = kirc_complete_ids))
kirc_com_snf_rand <-  as.data.frame(cbind(label = kirc_union_snf_rand, id = kirc_union_ids))
kirc_com_snf_self <- as.data.frame(cbind(label = kirc_com_snf_self, id = kirc_complete_ids))
kirc_com_snf_med <- as.data.frame(cbind(label = kirc_com_snf_med, id = kirc_complete_ids))
kirc_com_snf_reg <- as.data.frame(cbind(label = kirc_com_snf_reg, id = kirc_complete_ids))
kirc_union_snf_knn <- as.data.frame(cbind(label = kirc_union_snf_knn, id = kirc_union_ids))
kirc_union_snf_lls <- as.data.frame(cbind(label = kirc_union_snf_lls, id = kirc_union_ids))
kirc_union_snf_lsa <- as.data.frame(cbind(label = kirc_union_snf_lsa, id = kirc_union_ids))
kirc_union_snf_rand <-  as.data.frame(cbind(label = kirc_union_snf_rand, id = kirc_union_ids))
kirc_union_snf_self <-  as.data.frame(cbind(label = kirc_union_snf_self, id = kirc_union_ids))
kirc_union_snf_med <-  as.data.frame(cbind(label = kirc_union_snf_med, id = kirc_union_ids))
kirc_union_snf_reg <-  as.data.frame(cbind(label = kirc_union_snf_reg, id = kirc_union_ids))

findCommonCluster <- function(data1, data2) {
  numClus <- 5
  data_intersect_union <- data2[data2$id %in% data1$id,]
  data_intersect <- left_join(data_intersect_union, data1, by = 'id')
  common_clusters <- matrix(,0,5)
  cluster_compare <- matrix(,0,5)
  for(i in 1:numClus){
    sub_union <- data_intersect[data_intersect$label.x == i,]
    for(j in 1:numClus){
      sub_intersect <- data_intersect[data_intersect$label.y == j,]
      cluster_compare[j] <- round((sum(sub_union$id %in% sub_intersect$id)/nrow(sub_union))*100, 2)
    }
    common_clusters <- rbind(common_clusters, cluster_compare)
  }
  return(common_clusters)
}

kirc_hier_knn_common <- findCommonCluster(kirc_com_hier_knn, kirc_union_hier_knn)
kirc_hier_lls_common <- findCommonCluster(kirc_com_hier_lls, kirc_union_hier_lls)
kirc_hier_lsa_common <- findCommonCluster(kirc_com_hier_lsa, kirc_union_hier_lsa)
kirc_hier_rand_common <- findCommonCluster(kirc_com_hier_rand, kirc_union_hier_rand)

kirc_iclus_knn_common <- findCommonCluster(kirc_com_iclus_knn, kirc_union_iclus_knn)
kirc_iclus_lls_common <- findCommonCluster(kirc_com_iclus_lls, kirc_union_iclus_lls)
kirc_iclus_lsa_common <- findCommonCluster(kirc_com_iclus_lsa, kirc_union_iclus_lsa)
kirc_iclus_rand_common <- findCommonCluster(kirc_com_iclus_rand, kirc_union_iclus_rand)

kirc_snf_knn_common <- findCommonCluster(kirc_com_snf_knn, kirc_union_snf_knn)
kirc_snf_lls_common <- findCommonCluster(kirc_com_snf_lls, kirc_union_snf_lls)
kirc_snf_lsa_common <- findCommonCluster(kirc_com_snf_lsa, kirc_union_snf_lsa)
kirc_snf_rand_common <- findCommonCluster(kirc_com_snf_rand, kirc_union_snf_rand)
kirc_snf_self_common <- findCommonCluster(kirc_com_snf_self, kirc_union_snf_self)
kirc_snf_reg_common <- findCommonCluster(kirc_com_snf_reg, kirc_union_snf_reg)
kirc_snf_med_common <- findCommonCluster(kirc_com_snf_med, kirc_union_snf_med)
##########################################################################################################
# lihc, compare clusters from complete, compete imputed, and union

# get lihc IDs
lihc_complete_ids <- lihc_ids[[1]]
lihc_union_ids <- lihc_ids[[2]]

# Get labels for each clustering type, from the intersection for lihc
lihc_com_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_1.txt'))))
lihc_com_iclus <- as.factor(t(read.table(paste0(completeFolder, '/1_2.txt'))))
lihc_com_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_3.txt'))))

# get labels for lihc intersection imputed
# hier
lihc_com_hier_knn <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_knn_hier_comp.txt'))))
lihc_com_hier_lls <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_lls_hier_comp.txt'))))
lihc_com_hier_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_lsa_hier_comp.txt'))))
lihc_com_hier_rand <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_rand_hier_comp.txt'))))

#iclust
lihc_com_iclus_knn <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_knn_iclust_comp.txt'))))
lihc_com_iclus_lls <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_lls_iclust_comp.txt'))))
lihc_com_iclus_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_lsa_iclust_comp.txt'))))
lihc_com_iclus_rand <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_rand_iclust_comp.txt'))))

#SNF
lihc_com_snf_knn <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_knn_snf_comp.txt'))))
lihc_com_snf_lls <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_lls_snf_comp.txt'))))
lihc_com_snf_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_lsa_snf_comp.txt'))))
lihc_com_snf_rand <- as.factor(t(read.table(paste0(imputeFolder, '/lihc', '/lihc_rand_snf_comp.txt'))))
lihc_com_snf_self <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/lihc_self_com_labels.txt'))))
lihc_com_snf_med <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/lihc_med_com_labels.txt' ))))
lihc_com_snf_reg <- as.factor(t(read.table(paste0(similarityFolder, '/lihc', '/lihc_reg_com_labels.txt' ))))

###################################################################################################
# get labels for lihc union

# Hier
lihc_union_hier_knn  <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_knn_hier_union.txt'))))
lihc_union_hier_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_lls_hier_union.txt'))))
lihc_union_hier_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_lsa_hier_union.txt'))))
lihc_union_hier_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_rand_hier_union.txt'))))

#iclust
lihc_union_iclus_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_knn_iclust_union.txt'))))
lihc_union_iclus_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_lls_iclust_union.txt'))))
lihc_union_iclus_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_lsa_iclust_union.txt'))))
lihc_union_iclus_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_rand_iclust_union.txt'))))


# SNF
lihc_union_snf_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_knn_snf_union.txt'))))
lihc_union_snf_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_lls_snf_union.txt'))))
lihc_union_snf_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_lsa_snf_union.txt'))))
lihc_union_snf_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/lihc/', 'lihc_rand_snf_union.txt'))))
lihc_union_snf_self <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/lihc_self_union_labels.txt' ))))
lihc_union_snf_med <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/lihc_med_union_labels.txt' ))))
lihc_union_snf_reg <- as.factor(t(read.table(paste0(similarityOrigFolder, '/lihc', '/lihc_reg_union_labels.txt' ))))

##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical 
lihc_hier_stat <- rbind(
  append('lihc_com_hier_knn', summary(as.factor(lihc_com_hier_knn))),
  append('lihc_com_hier_lls', summary(as.factor(lihc_com_hier_lls))),
  append('lihc_com_hier_lsa', summary(as.factor(lihc_com_hier_lsa))),
  append('lihc_com_hier_rand', summary(as.factor(lihc_com_hier_rand))),
  append('lihc_union_hier_knn', summary(as.factor(lihc_union_hier_knn))),
  append('lihc_union_hier_lls', summary(as.factor(lihc_union_hier_lls))),
  append('lihc_union_hier_lsa', summary(as.factor(lihc_union_hier_lsa))),
  append('lihc_union_hier_rand', summary(as.factor(lihc_union_hier_rand)))
  
)

# Icluster
lihc_iclust_stat <- rbind(
  append('lihc_com_iclust_knn', summary(as.factor(lihc_com_iclus_knn))),
  append('lihc_com_iclust_lls', summary(as.factor(lihc_com_iclus_lls))),
  append('lihc_com_iclust_lsa', summary(as.factor(lihc_com_iclus_lsa))),
  append('lihc_com_iclust_rand', summary(as.factor(lihc_com_iclus_rand))),
  append('lihc_union_iclust_knn', summary(as.factor(lihc_union_iclus_knn))),
  append('lihc_union_iclust_lls', summary(as.factor(lihc_union_iclus_lls))),
  append('lihc_union_iclust_lsa', summary(as.factor(lihc_union_iclus_lsa))),
  append('lihc_union_iclust_rand', summary(as.factor(lihc_union_iclus_rand)))
  
)


# SNF
lihc_snf_stat <- rbind(
  append('lihc_com_snf', summary(as.factor(lihc_com_snf))),
  append('lihc_com_snf_self', summary(as.factor(lihc_com_snf_self))),
  append('lihc_com_snf_med', summary(as.factor(lihc_com_snf_med))),
  append('lihc_com_snf_reg', summary(as.factor(lihc_com_snf_reg))),
  append('lihc_union_snf_self', summary(as.factor(lihc_union_snf_self))),
  append('lihc_union_snf_med', summary(as.factor(lihc_union_snf_med))),
  append('lihc_union_snf_reg', summary(as.factor(lihc_union_snf_reg)))
)



# Merge ids with labels and get clinical information for each cluster

# look at groups of ids in complete, complete_impute, and union
# join ids and then join each label groupd
lihc_com_hier <- as.data.frame(cbind(label = lihc_com_hier, id = lihc_complete_ids))
lihc_com_hier_knn <- as.data.frame(cbind(label = lihc_com_hier_knn, id = lihc_complete_ids))
lihc_com_hier_lls <- as.data.frame(cbind(label = lihc_com_hier_lls, id = lihc_complete_ids))
lihc_com_hier_lsa <- as.data.frame(cbind(label = lihc_com_hier_lsa, id = lihc_complete_ids))
lihc_com_hier_rand <-  as.data.frame(cbind(label = lihc_union_hier_rand, id = lihc_union_ids))
lihc_union_hier_knn <- as.data.frame(cbind(label = lihc_union_hier_knn, id = lihc_union_ids))
lihc_union_hier_lls <- as.data.frame(cbind(label = lihc_union_hier_lls, id = lihc_union_ids))
lihc_union_hier_lsa <- as.data.frame(cbind(label = lihc_union_hier_lsa, id = lihc_union_ids))
lihc_union_hier_rand <-  as.data.frame(cbind(label = lihc_union_hier_rand, id = lihc_union_ids))

lihc_com_iclus <- as.data.frame(cbind(label = lihc_com_iclus, id = lihc_complete_ids))
lihc_com_iclus_knn <- as.data.frame(cbind(label = lihc_com_iclus_knn, id = lihc_complete_ids))
lihc_com_iclus_lls <- as.data.frame(cbind(label = lihc_com_iclus_lls, id = lihc_complete_ids))
lihc_com_iclus_lsa <- as.data.frame(cbind(label = lihc_com_iclus_lsa, id = lihc_complete_ids))
lihc_com_iclus_rand <-  as.data.frame(cbind(label = lihc_union_iclus_rand, id = lihc_union_ids))
lihc_union_iclus_knn <- as.data.frame(cbind(label = lihc_union_iclus_knn, id = lihc_union_ids))
lihc_union_iclus_lls <- as.data.frame(cbind(label = lihc_union_iclus_lls, id = lihc_union_ids))
lihc_union_iclus_lsa <- as.data.frame(cbind(label = lihc_union_iclus_lsa, id = lihc_union_ids))
lihc_union_iclus_rand <-  as.data.frame(cbind(label = lihc_union_iclus_rand, id = lihc_union_ids))

lihc_com_snf <- as.data.frame(cbind(label = lihc_com_snf, id = lihc_complete_ids))
lihc_com_snf_knn <- as.data.frame(cbind(label = lihc_com_snf_knn, id = lihc_complete_ids))
lihc_com_snf_lls <- as.data.frame(cbind(label = lihc_com_snf_lls, id = lihc_complete_ids))
lihc_com_snf_lsa <- as.data.frame(cbind(label = lihc_com_snf_lsa, id = lihc_complete_ids))
lihc_com_snf_rand <-  as.data.frame(cbind(label = lihc_union_snf_rand, id = lihc_union_ids))
lihc_com_snf_self <- as.data.frame(cbind(label = lihc_com_snf_self, id = lihc_complete_ids))
lihc_com_snf_med <- as.data.frame(cbind(label = lihc_com_snf_med, id = lihc_complete_ids))
lihc_com_snf_reg <- as.data.frame(cbind(label = lihc_com_snf_reg, id = lihc_complete_ids))
lihc_union_snf_knn <- as.data.frame(cbind(label = lihc_union_snf_knn, id = lihc_union_ids))
lihc_union_snf_lls <- as.data.frame(cbind(label = lihc_union_snf_lls, id = lihc_union_ids))
lihc_union_snf_lsa <- as.data.frame(cbind(label = lihc_union_snf_lsa, id = lihc_union_ids))
lihc_union_snf_rand <-  as.data.frame(cbind(label = lihc_union_snf_rand, id = lihc_union_ids))
lihc_union_snf_self <-  as.data.frame(cbind(label = lihc_union_snf_self, id = lihc_union_ids))
lihc_union_snf_med <-  as.data.frame(cbind(label = lihc_union_snf_med, id = lihc_union_ids))
lihc_union_snf_reg <-  as.data.frame(cbind(label = lihc_union_snf_reg, id = lihc_union_ids))

findCommonCluster <- function(data1, data2) {
  numClus <- 5
  data_intersect_union <- data2[data2$id %in% data1$id,]
  data_intersect <- left_join(data_intersect_union, data1, by = 'id')
  common_clusters <- matrix(,0,5)
  cluster_compare <- matrix(,0,5)
  for(i in 1:numClus){
    sub_union <- data_intersect[data_intersect$label.x == i,]
    for(j in 1:numClus){
      sub_intersect <- data_intersect[data_intersect$label.y == j,]
      cluster_compare[j] <- round((sum(sub_union$id %in% sub_intersect$id)/nrow(sub_union))*100, 2)
    }
    common_clusters <- rbind(common_clusters, cluster_compare)
  }
  return(common_clusters)
}

lihc_hier_knn_common <- findCommonCluster(lihc_com_hier_knn, lihc_union_hier_knn)
lihc_hier_lls_common <- findCommonCluster(lihc_com_hier_lls, lihc_union_hier_lls)
lihc_hier_lsa_common <- findCommonCluster(lihc_com_hier_lsa, lihc_union_hier_lsa)
lihc_hier_rand_common <- findCommonCluster(lihc_com_hier_rand, lihc_union_hier_rand)

lihc_iclus_knn_common <- findCommonCluster(lihc_com_iclus_knn, lihc_union_iclus_knn)
lihc_iclus_lls_common <- findCommonCluster(lihc_com_iclus_lls, lihc_union_iclus_lls)
lihc_iclus_lsa_common <- findCommonCluster(lihc_com_iclus_lsa, lihc_union_iclus_lsa)
lihc_iclus_rand_common <- findCommonCluster(lihc_com_iclus_rand, lihc_union_iclus_rand)

lihc_snf_knn_common <- findCommonCluster(lihc_com_snf_knn, lihc_union_snf_knn)
lihc_snf_lls_common <- findCommonCluster(lihc_com_snf_lls, lihc_union_snf_lls)
lihc_snf_lsa_common <- findCommonCluster(lihc_com_snf_lsa, lihc_union_snf_lsa)
lihc_snf_rand_common <- findCommonCluster(lihc_com_snf_rand, lihc_union_snf_rand)
lihc_snf_self_common <- findCommonCluster(lihc_com_snf_self, lihc_union_snf_self)
lihc_snf_reg_common <- findCommonCluster(lihc_com_snf_reg, lihc_union_snf_reg)
lihc_snf_med_common <- findCommonCluster(lihc_com_snf_med, lihc_union_snf_med)
##########################################################################################################
# luad, compare clusters from complete, compete imputed, and union

# get luad IDs
luad_complete_ids <- luad_ids[[1]]
luad_union_ids <- luad_ids[[2]]

# Get labels for each clustering type, from the intersection for luad
luad_com_hier <- as.factor(t(read.table(paste0(completeFolder, '/1_1.txt'))))
luad_com_iclus <- as.factor(t(read.table(paste0(completeFolder, '/1_2.txt'))))
luad_com_snf <- as.factor(t(read.table(paste0(completeFolder, '/1_3.txt'))))

# get labels for luad intersection imputed
# hier
luad_com_hier_knn <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_knn_hier_comp.txt'))))
luad_com_hier_lls <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_lls_hier_comp.txt'))))
luad_com_hier_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_lsa_hier_comp.txt'))))
luad_com_hier_rand <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_rand_hier_comp.txt'))))

#iclust
luad_com_iclus_knn <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_knn_iclust_comp.txt'))))
luad_com_iclus_lls <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_lls_iclust_comp.txt'))))
luad_com_iclus_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_lsa_iclust_comp.txt'))))
luad_com_iclus_rand <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_rand_iclust_comp.txt'))))

#SNF
luad_com_snf_knn <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_knn_snf_comp.txt'))))
luad_com_snf_lls <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_lls_snf_comp.txt'))))
luad_com_snf_lsa <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_lsa_snf_comp.txt'))))
luad_com_snf_rand <- as.factor(t(read.table(paste0(imputeFolder, '/luad', '/luad_rand_snf_comp.txt'))))
luad_com_snf_self <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/luad_self_com_labels.txt'))))
luad_com_snf_med <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/luad_med_com_labels.txt' ))))
luad_com_snf_reg <- as.factor(t(read.table(paste0(similarityFolder, '/luad', '/luad_reg_com_labels.txt' ))))

###################################################################################################
# get labels for luad union

# Hier
luad_union_hier_knn  <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_knn_hier_union.txt'))))
luad_union_hier_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_lls_hier_union.txt'))))
luad_union_hier_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_lsa_hier_union.txt'))))
luad_union_hier_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_rand_hier_union.txt'))))

#iclust
luad_union_iclus_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_knn_iclust_union.txt'))))
luad_union_iclus_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_lls_iclust_union.txt'))))
luad_union_iclus_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_lsa_iclust_union.txt'))))
luad_union_iclus_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_rand_iclust_union.txt'))))


# SNF
luad_union_snf_knn <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_knn_snf_union.txt'))))
luad_union_snf_lls <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_lls_snf_union.txt'))))
luad_union_snf_lsa <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_lsa_snf_union.txt'))))
luad_union_snf_rand <- as.factor(t(read.table(paste0(imputeOrigFolder, '/luad/', 'luad_rand_snf_union.txt'))))
luad_union_snf_self <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/luad_self_union_labels.txt' ))))
luad_union_snf_med <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/luad_med_union_labels.txt' ))))
luad_union_snf_reg <- as.factor(t(read.table(paste0(similarityOrigFolder, '/luad', '/luad_reg_union_labels.txt' ))))

##########################################################################################################
# Compare labels for each clustering type 

# Hierarchical 
luad_hier_stat <- rbind(
  append('luad_com_hier_knn', summary(as.factor(luad_com_hier_knn))),
  append('luad_com_hier_lls', summary(as.factor(luad_com_hier_lls))),
  append('luad_com_hier_lsa', summary(as.factor(luad_com_hier_lsa))),
  append('luad_com_hier_rand', summary(as.factor(luad_com_hier_rand))),
  append('luad_union_hier_knn', summary(as.factor(luad_union_hier_knn))),
  append('luad_union_hier_lls', summary(as.factor(luad_union_hier_lls))),
  append('luad_union_hier_lsa', summary(as.factor(luad_union_hier_lsa))),
  append('luad_union_hier_rand', summary(as.factor(luad_union_hier_rand)))
  
)

# Icluster
luad_iclust_stat <- rbind(
  append('luad_com_iclust_knn', summary(as.factor(luad_com_iclus_knn))),
  append('luad_com_iclust_lls', summary(as.factor(luad_com_iclus_lls))),
  append('luad_com_iclust_lsa', summary(as.factor(luad_com_iclus_lsa))),
  append('luad_com_iclust_rand', summary(as.factor(luad_com_iclus_rand))),
  append('luad_union_iclust_knn', summary(as.factor(luad_union_iclus_knn))),
  append('luad_union_iclust_lls', summary(as.factor(luad_union_iclus_lls))),
  append('luad_union_iclust_lsa', summary(as.factor(luad_union_iclus_lsa))),
  append('luad_union_iclust_rand', summary(as.factor(luad_union_iclus_rand)))
  
)


# SNF
luad_snf_stat <- rbind(
  append('luad_com_snf', summary(as.factor(luad_com_snf))),
  append('luad_com_snf_self', summary(as.factor(luad_com_snf_self))),
  append('luad_com_snf_med', summary(as.factor(luad_com_snf_med))),
  append('luad_com_snf_reg', summary(as.factor(luad_com_snf_reg))),
  append('luad_union_snf_self', summary(as.factor(luad_union_snf_self))),
  append('luad_union_snf_med', summary(as.factor(luad_union_snf_med))),
  append('luad_union_snf_reg', summary(as.factor(luad_union_snf_reg)))
)



# Merge ids with labels and get clinical information for each cluster

# look at groups of ids in complete, complete_impute, and union
# join ids and then join each label groupd
luad_com_hier <- as.data.frame(cbind(label = luad_com_hier, id = luad_complete_ids))
luad_com_hier_knn <- as.data.frame(cbind(label = luad_com_hier_knn, id = luad_complete_ids))
luad_com_hier_lls <- as.data.frame(cbind(label = luad_com_hier_lls, id = luad_complete_ids))
luad_com_hier_lsa <- as.data.frame(cbind(label = luad_com_hier_lsa, id = luad_complete_ids))
luad_com_hier_rand <-  as.data.frame(cbind(label = luad_union_hier_rand, id = luad_union_ids))
luad_union_hier_knn <- as.data.frame(cbind(label = luad_union_hier_knn, id = luad_union_ids))
luad_union_hier_lls <- as.data.frame(cbind(label = luad_union_hier_lls, id = luad_union_ids))
luad_union_hier_lsa <- as.data.frame(cbind(label = luad_union_hier_lsa, id = luad_union_ids))
luad_union_hier_rand <-  as.data.frame(cbind(label = luad_union_hier_rand, id = luad_union_ids))

luad_com_iclus <- as.data.frame(cbind(label = luad_com_iclus, id = luad_complete_ids))
luad_com_iclus_knn <- as.data.frame(cbind(label = luad_com_iclus_knn, id = luad_complete_ids))
luad_com_iclus_lls <- as.data.frame(cbind(label = luad_com_iclus_lls, id = luad_complete_ids))
luad_com_iclus_lsa <- as.data.frame(cbind(label = luad_com_iclus_lsa, id = luad_complete_ids))
luad_com_iclus_rand <-  as.data.frame(cbind(label = luad_union_iclus_rand, id = luad_union_ids))
luad_union_iclus_knn <- as.data.frame(cbind(label = luad_union_iclus_knn, id = luad_union_ids))
luad_union_iclus_lls <- as.data.frame(cbind(label = luad_union_iclus_lls, id = luad_union_ids))
luad_union_iclus_lsa <- as.data.frame(cbind(label = luad_union_iclus_lsa, id = luad_union_ids))
luad_union_iclus_rand <-  as.data.frame(cbind(label = luad_union_iclus_rand, id = luad_union_ids))

luad_com_snf <- as.data.frame(cbind(label = luad_com_snf, id = luad_complete_ids))
luad_com_snf_knn <- as.data.frame(cbind(label = luad_com_snf_knn, id = luad_complete_ids))
luad_com_snf_lls <- as.data.frame(cbind(label = luad_com_snf_lls, id = luad_complete_ids))
luad_com_snf_lsa <- as.data.frame(cbind(label = luad_com_snf_lsa, id = luad_complete_ids))
luad_com_snf_rand <-  as.data.frame(cbind(label = luad_union_snf_rand, id = luad_union_ids))
luad_com_snf_self <- as.data.frame(cbind(label = luad_com_snf_self, id = luad_complete_ids))
luad_com_snf_med <- as.data.frame(cbind(label = luad_com_snf_med, id = luad_complete_ids))
luad_com_snf_reg <- as.data.frame(cbind(label = luad_com_snf_reg, id = luad_complete_ids))
luad_union_snf_knn <- as.data.frame(cbind(label = luad_union_snf_knn, id = luad_union_ids))
luad_union_snf_lls <- as.data.frame(cbind(label = luad_union_snf_lls, id = luad_union_ids))
luad_union_snf_lsa <- as.data.frame(cbind(label = luad_union_snf_lsa, id = luad_union_ids))
luad_union_snf_rand <-  as.data.frame(cbind(label = luad_union_snf_rand, id = luad_union_ids))
luad_union_snf_self <-  as.data.frame(cbind(label = luad_union_snf_self, id = luad_union_ids))
luad_union_snf_med <-  as.data.frame(cbind(label = luad_union_snf_med, id = luad_union_ids))
luad_union_snf_reg <-  as.data.frame(cbind(label = luad_union_snf_reg, id = luad_union_ids))

findCommonCluster <- function(data1, data2) {
  numClus <- 5
  data_intersect_union <- data2[data2$id %in% data1$id,]
  data_intersect <- left_join(data_intersect_union, data1, by = 'id')
  common_clusters <- matrix(,0,5)
  cluster_compare <- matrix(,0,5)
  for(i in 1:numClus){
    sub_union <- data_intersect[data_intersect$label.x == i,]
    for(j in 1:numClus){
      sub_intersect <- data_intersect[data_intersect$label.y == j,]
      cluster_compare[j] <- round((sum(sub_union$id %in% sub_intersect$id)/nrow(sub_union))*100, 2)
    }
    common_clusters <- rbind(common_clusters, cluster_compare)
  }
  return(common_clusters)
}

luad_hier_knn_common <- findCommonCluster(luad_com_hier_knn, luad_union_hier_knn)
luad_hier_lls_common <- findCommonCluster(luad_com_hier_lls, luad_union_hier_lls)
luad_hier_lsa_common <- findCommonCluster(luad_com_hier_lsa, luad_union_hier_lsa)
luad_hier_rand_common <- findCommonCluster(luad_com_hier_rand, luad_union_hier_rand)

luad_iclus_knn_common <- findCommonCluster(luad_com_iclus_knn, luad_union_iclus_knn)
luad_iclus_lls_common <- findCommonCluster(luad_com_iclus_lls, luad_union_iclus_lls)
luad_iclus_lsa_common <- findCommonCluster(luad_com_iclus_lsa, luad_union_iclus_lsa)
luad_iclus_rand_common <- findCommonCluster(luad_com_iclus_rand, luad_union_iclus_rand)

luad_snf_knn_common <- findCommonCluster(luad_com_snf_knn, luad_union_snf_knn)
luad_snf_lls_common <- findCommonCluster(luad_com_snf_lls, luad_union_snf_lls)
luad_snf_lsa_common <- findCommonCluster(luad_com_snf_lsa, luad_union_snf_lsa)
luad_snf_rand_common <- findCommonCluster(luad_com_snf_rand, luad_union_snf_rand)
luad_snf_self_common <- findCommonCluster(luad_com_snf_self, luad_union_snf_self)
luad_snf_reg_common <- findCommonCluster(luad_com_snf_reg, luad_union_snf_reg)
luad_snf_med_common <- findCommonCluster(luad_com_snf_med, luad_union_snf_med)

