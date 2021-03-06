##### For evaluating and ranking survival metrics for original imputation and similarity

# Scripts/06_Two_Thousand_Features/evaluate_imputation/Results

# Evalute the results of the evaluate_imputation tests

######################################################################
# Initialize folders
home_folder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
project_folder <- paste(home_folder, 'Projects/SNF/NM_2015', sep = '/')
scripts_folder <- paste(project_folder, "Scripts", "06_Combat", sep = "/")
plots_folder <- paste(scripts_folder, "analyze_results/Plots", sep ="/")
results_folder <- paste(project_folder, 'Scripts/06_Results', sep = '/')


######################################################################
# Initialize fixed variables
cancerTypes <- 'LUSC'
dataTypes <- c("methyl", "mirna", "mrna")
imputeTypes <- c("KNN", "LLS", "LSA", "Rand")
clusterTypes <- c("Hierarchical", "iCluster", "SNF")
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
clusterCols <- c("impute", "seed", "cluster", "set", "acc",
                 "nmi", "pval", "ci", "se.std", "pvalcox", "coefcox", "con_index_ci", "con_index_p", 
                 "bias_corrected_c_index")

# cluster <- loadResults("evaluate_imputation", "clustering.txt",
#                        clusterCols)

originalCluster <- loadResults("evaluate_original_imputation",
                               "clustering.txt", clusterCols[-2])

# Load the similarity results
similarityCols <- c("seed", "similarity", "set", "acc",
                    "nmi", "pval", "ci", "se.std", "pvalcox", "coefcox", "con_index_ci", "con_index_p", 
                    "bias_corrected_c_index")

# similarity <- loadResults("evaluate_similarity", "similarity.txt",
#                           similarityCols)
originalSimilarity <- loadResults("evaluate_original_similarity",
                                  "similarity.txt",
                                  similarityCols[-1])

######################################################################
# Transform pval to -log10(pval)
# 
# cluster$pval <- -log10(cluster$pval)
# cluster$pvalcox <- -log10(cluster$pvalcox)
# cluster$con_index_p <- -log10(cluster$con_index_p)
originalCluster$pval <- -log10(originalCluster$pval)
originalCluster$pvalcox <- -log10(originalCluster$pvalcox)
originalCluster$con_index_p <- -log10(originalCluster$con_index_p)
# similarity$pval <- -log10(similarity$pval)
# similarity$pvalcox <- -log10(similarity$pvalcox)
# similarity$con_index_p <- -log10(similarity$con_index_p)
originalSimilarity$pval <- -log10(originalSimilarity$pval)
originalSimilarity$pvalcox <- -log10(originalSimilarity$pvalcox)
originalSimilarity$con_index_p <- -log10(originalSimilarity$con_index_p)

############################################################################
# Rank original clustering and similarity resutls

# Combine the similarity and clustering data sets
# Similarity:
# -add "cluster" column, with "SNF" in each entry
# -rename "similarity" column to "impute"
newSimilarity <- originalSimilarity
colnames(newSimilarity)[1] <- "impute"
newSimilarity$cluster <- "SNF"
newSimilarity$impute <- similarityTypes[newSimilarity$impute]

newCluster <- originalCluster
newCluster$cluster <- clusterTypes[newCluster$cluster]
newCluster$impute <- imputeTypes[newCluster$impute]

combinedScores <- rbind(newCluster, newSimilarity)
combinedScores <- subset(combinedScores, set==1,
                         c(impute, cluster, acc, nmi,
                           pval, ci, se.std, pvalcox, coefcox, con_index_ci, con_index_p, 
                           bias_corrected_c_index))

combinedScores$method <- interaction(combinedScores$cluster,
                                     combinedScores$impute, drop=TRUE)
methodTypes <- levels(combinedScores$method)
combinedScores$method <- as.numeric(combinedScores$method)

# Remove Inf pval's
combinedScores$pval[is.infinite(combinedScores$pval)] <- NA

# Write combinedScores to results folder 
write.csv(combinedScores, paste0(results_folder, '/scoresCombatOrigDup.csv'))


# Compare if the mean of method 1 is significantly greater than the
# mean of method 2, for every pair of methods and cancer type


testScores <- matrix(, 0, 5)


for (met1 in 1:length(methodTypes)) {
  score1 <- subset(combinedScores, method==met1)
  accVector <- rep.int(0, length(methodTypes))
  nmiVector <- accVector
  pvalVector <- accVector
  ciVector <- accVector
  for (met2 in (1:length(methodTypes))[-met1]) {
    score2 <- subset(combinedScores, method==met2)
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

write.table(testScores, "testScoresOriginal.txt", row.names=FALSE)
