##### For evaluating and ranking survival metrics for original imputation and similarity

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
scripts_folder <- paste(project_folder, "Scripts", "06_Int", sep = "/")
plots_folder <- paste(scripts_folder, "analyze_results/Plots", sep ="/")
results_folder <- paste(project_folder, 'Scripts/06_Results', sep = '/')


######################################################################
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD") # Removed LUSC
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
clusterCols <- c("cancer", "impute", "seed", "cluster", "set", "acc",
                 "nmi", "pval", "ci", "intAcc", "intNmi", "intPval",
                 "intCi", "runTime")
imputeCols <- c("cancer", "impute", "seed", "methyl", "mirna", "mrna",
                "runtime")

originalCluster <- loadResults("evaluate_original_imputation",
                               "clustering.txt", clusterCols[-3])

# Load the similarity results
similarityCols <- c("cancer", "seed", "similarity", "set", "acc",
                    "nmi", "pval", "ci", "intAcc", "intNmi",
                    "intPval", "intCi", "runTime")

originalSimilarity <- loadResults("evaluate_original_similarity",
                                  "similarity.txt",
                                  similarityCols[-2])

######################################################################
# Transform pval to -log10(pval)


originalCluster$actual_pval <- -log10(originalCluster$pval)
originalCluster$pval <- -log10(originalCluster$pval)

originalSimilarity$actual_pval <- originalSimilarity$pval
originalSimilarity$pval <- -log10(originalSimilarity$pval)
############################################################################
# Rank original clustering and similarity resutls

# Combine the similarity and clustering data sets
# Similarity:
# -add "cluster" column, with "SNF" in each entry
# -rename "similarity" column to "impute"
newSimilarity <- originalSimilarity
colnames(newSimilarity)[similarityCols[-2]=="similarity"] <- "impute"
newSimilarity$cluster <- "SNF"
newSimilarity$impute <- similarityTypes[newSimilarity$impute]

newCluster <- originalCluster
newCluster$cluster <- clusterTypes[newCluster$cluster]
newCluster$impute <- imputeTypes[newCluster$impute]

combinedScores <- rbind(newCluster, newSimilarity)
combinedScores <- subset(combinedScores, set==1,
                         c(cancer, impute, cluster, intAcc, intNmi,
                           intPval, intCi))
combinedScores$method <- interaction(combinedScores$cluster,
                                     combinedScores$impute, drop=TRUE)
methodTypes <- levels(combinedScores$method)
combinedScores$method <- as.numeric(combinedScores$method)

# Remove Inf pval's
combinedScores$intPval[is.infinite(combinedScores$intPval)] <- NA

# Write combinedScores to results folder 
write.csv(combinedScores, paste0(results_folder, '/scoresTwoThousandOrigInt.csv'))

# Compare if the mean of method 1 is significantly greater than the
# mean of method 2, for every pair of methods and cancer type
xGreaterY <- function(x, y) {
  t.test(x, y, alternative="greater", paired=TRUE, na.action = na.omit)$p.value < 0.05
}

testScores <- matrix(, 0, 6)

for (canc in 1:length(cancerTypes)) {
  for (met1 in 1:length(methodTypes)) {
    score1 <- subset(combinedScores, cancer==canc&method==met1)
    intAccVector <- rep.int(0, length(methodTypes))
    intNmiVector <- intAccVector
    intPvalVector <- intAccVector
    intCiVector <- intAccVector
    for (met2 in (1:length(methodTypes))[-met1]) {
      score2 <- subset(combinedScores, cancer==canc&method==met2)
      intAccVector[met2] <- score1$intAcc > score2$intAcc
      intNmiVector[met2] <- score1$intNmi > score2$intNmi
      intPvalVector[met2] <-score1$intPval > score2$intPval
      intCiVector[met2] <-  score1$intCi > score2$intCi
    }
    testScores <- rbind(testScores,
                        c(canc, met1, sum(intAccVector), sum(intNmiVector),
                          sum(intPvalVector), sum(intCiVector)))
  }
}

colnames(testScores) <- c("cancer", "method", "intAcc", "intNmi", "intPval",
                          "intCi")
testScores <- as.data.frame(testScores)
testScores$cancer <- cancerTypes[testScores$cancer]
testScores$method <- methodTypes[testScores$method]

write.table(testScores, "testScoresOriginal.txt", row.names=FALSE)
