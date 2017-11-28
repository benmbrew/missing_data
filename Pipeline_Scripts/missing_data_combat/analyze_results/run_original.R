# Need to evaluate the clustering results

# Scripts/06_Two_Thousand_Features/evaluate_imputation/Results

# Evalute the results of the evaluate_imputation tests

######################################################################
# Load libraries
# library(dplyr)
# library(reshape2)

######################################################################
# Initialize folders
home_folder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
project_folder <- paste(home_folder, 'Projects/SNF/NM_2015', sep = '/')
scripts_folder <- paste(project_folder, "Scripts", "missing_data_combat", sep = "/")
plots_folder <- paste(scripts_folder, "analyze_results/Plots", sep ="/")
results_folder <- paste(project_folder, 'Scripts/06_Results', sep = '/')


######################################################################
# Initialize fixed variables
dataTypes <- c("methyl", "mirna", "mrna")
imputeTypes <- c("KNN", "LLS", "LSA", "Rand")
setTypes <- c("Union", "Intersect")
similarityTypes <- c("Self", "Median", "Regres")

######################################################################
# Load the results
# Accuracy and NMI metrics are meaningless in the original results
# because they were compared to randomly generated labels
loadResults <- function(test_folder, file_name, colNames=NULL) {
  file_path <- paste(scripts_folder, test_folder, "Results", file_name,
                     sep="/")
  read.table(file_path, header=FALSE, col.names=colNames)
}

# Load the imputation results
clusterCols <- c("impute", "cluster_size", "acc",
                 "nmi", "pval", "ci", "se.std", "pvalcox", "coefcox", "con_index_ci", "con_index_p", 
                 "bias_corrected_c_index")

cluster <- loadResults("evaluate_original_imputation", "clustering.txt",
                       clusterCols)

# Load the similarity results
similarityCols <- c("similarity", "cluster_size", "acc",
                    "nmi", "pval", "ci", "se.std", "pvalcox", "coefcox", "con_index_ci",
                    "con_index_p", "bias_corrected_c_index")


similarity <- loadResults("evaluate_original_similarity", "similarity.txt",
                          similarityCols)

######################################################################
# Transform pval to -log10(pval)


cluster$pval <- -log10(cluster$pval)
cluster$pvalcox <- -log10(cluster$pvalcox)
cluster$con_index_p <- -log10(cluster$con_index_p)

similarity$pval <- -log10(similarity$pval)
similarity$pvalcox <- -log10(similarity$pvalcox)
similarity$con_index_p <- -log10(similarity$con_index_p)


# Combine the similarity and clustering data sets
# Similarity:
# -add "cluster" column, with "SNF" in each entry
# -rename "similarity" column to "impute"
newSimilarity <- similarity
colnames(newSimilarity)[similarityCols=="similarity"] <- "impute"
newSimilarity$cluster <- "SNF"
newSimilarity$impute <- similarityTypes[newSimilarity$impute]

newCluster <- cluster
newCluster$cluster <- "SNF"
newCluster$impute <- imputeTypes[newCluster$impute]

combinedScores <- rbind(newCluster, newSimilarity)
combinedScores <- combinedScores[,c('impute', 'cluster', 'cluster_size', 'acc', 'nmi', 
                                    'pval', 'ci', 'se.std', 'pvalcox')]
combinedScores$method <- interaction(combinedScores$cluster,
                                     combinedScores$impute, drop=TRUE)
methodTypes <- levels(combinedScores$method)
combinedScores$method <- as.numeric(combinedScores$method)

# Remove Inf pval's
combinedScores$pval[is.infinite(combinedScores$pval)] <- NA

# Write combinedScores to results folder 
write.csv(combinedScores, paste0(results_folder, '/missing_data_combat_union.csv'))


# Compare if the mean of method 1 is significantly greater than the
# mean of method 2, for every pair of methods and cancer type


byCluster <- function(data, cluster_size){
  
  testScores <- matrix(, 0, 5)
  
  data <- data[data$cluster_size == cluster_size,]
  
  for (met1 in 1:length(methodTypes)) {
    score1 <- subset(data, method==met1)
    accVector <- rep.int(0, length(methodTypes))
    nmiVector <- accVector
    pvalVector <- accVector
    ciVector <- accVector
    for (met2 in (1:length(methodTypes))[-met1]) {
      score2 <- subset(data, method==met2)
      accVector[met2] <- score1$acc > score2$acc
      nmiVector[met2] <- score1$nmi > score2$nmi
      pvalVector[met2] <-score1$pval > score2$pval
      ciVector[met2] <-  score1$ci > score2$ci
    }
    testScores <- rbind(testScores,
                        c(met1, sum(accVector), sum(nmiVector),
                          sum(pvalVector), sum(ciVector)))
  }
  
  
  colnames(testScores) <- c("method", "acc", "nmi", "pval",
                            "ci")
  testScores <- as.data.frame(testScores)
  testScores$method <- methodTypes[testScores$method]
  return(testScores)
}

combat5 <- byCluster(combinedScores, cluster_size= 1)
combat4 <- byCluster(combinedScores, cluster_size= 2)
combat3 <- byCluster(combinedScores, cluster_size= 3)
combat2 <- byCluster(combinedScores, cluster_size= 4)



write.table(testScores, "testScores.txt", row.names=FALSE)
