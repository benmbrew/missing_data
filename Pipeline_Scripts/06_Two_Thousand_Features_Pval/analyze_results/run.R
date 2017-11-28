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
scripts_folder <- paste(project_folder, "Scripts", "06_Two_Thousand_Features_Pval", sep = "/")
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
                 "nmi", "pval", "ci", "se.std", "pvalcox", "coefcox", "con_index_ci", "con_index_p", 
                 "bias_corrected_c_index", "intAcc", "intNmi", "intPval",
                 "intCi", "intSe.std", "intPvalx", "intCoefcox", "int_con_index_ci", 
                 "int_con_index_p", "int_bias_corrected_c_index", "runTime")

cluster <- loadResults("evaluate_imputation", "clustering.txt",
                       clusterCols)

originalCluster <- loadResults("evaluate_original_imputation",
                               "clustering.txt", clusterCols[-3])

# Load the similarity results
similarityCols <- c("cancer", "seed", "similarity", "set", "acc",
                    "nmi", "pval", "ci", "se.std", "pvalcox", "coefcox", "con_index_ci", "con_index_p", 
                    "bias_corrected_c_index", "intAcc", "intNmi",
                    "intPval", "intCi", "intSe.std", "intPvalx", "intCoefcox", "int_con_index_ci", 
                    "int_con_index_p", "int_bias_corrected_c_index", "runTime")

similarity <- loadResults("evaluate_similarity", "Similarity.txt",
                          similarityCols)
originalSimilarity <- loadResults("evaluate_original_similarity",
                                  "similarity.txt",
                                  similarityCols[-2])

######################################################################
# Transform pval to -log10(pval)

cluster$pval <- -log10(cluster$pval)
cluster$pvalcox <- -log10(cluster$pvalcox)
cluster$con_index_p <- -log10(cluster$con_index_p)
originalCluster$pval <- -log10(originalCluster$pval)
originalCluster$pvalcox <- -log10(originalCluster$pvalcox)
originalCluster$con_index_p <- -log10(originalCluster$con_index_p)
similarity$pval <- -log10(similarity$pval)
similarity$pvalcox <- -log10(similarity$pvalcox)
similarity$con_index_p <- -log10(similarity$con_index_p)
originalSimilarity$pval <- -log10(originalSimilarity$pval)
originalSimilarity$pvalcox <- -log10(originalSimilarity$pvalcox)
originalSimilarity$con_index_p <- -log10(originalSimilarity$con_index_p)


######################################################################
# Load helper functions

setPlotLayout <- function(xDim, yDim, byRow) {
  dims <- c(xDim, yDim)
  if (byRow) {
    par(mfrow=dims)
  } else {
    par(mfcol=dims)
  }
  par(cex=0.6)
  par(mar=c(0, 0, 0, 0), oma=c(4, 4, 2, 2))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
}

labelPlot <- function(title, xLabel, yLabel) {
  labels <- c(xLabel, yLabel, title)
  for (i in 1:length(labels)) {
    mtext(labels[i], side=i, outer=TRUE, cex=0.7,
          line=ifelse(i==3, 1, 2.2), col="grey20")
  }
}

######################################################################
# Plot the results

# 1. Plot each imputation method's error on all of the data types
# -Three box plots (one plot for each data type)
# -Three bars in each graph (one for each imputation method)
imputeBoxplot <- function(formula, canc, ylim, ny=10, label) {
  boxplot(formula=formula, data=impute, subset=cancer==canc,
          ylim=ylim, names=imputeTypes)
  grid(nx=NA, ny=ny, col="grey40")
  mtext(paste(cancerTypes[canc], label, sep=", "), side=3, line=-1,
        adj=0.1, cex=0.6, col='grey40')
}

pdf(file=paste0(plotsFolder, "/imputation.pdf"))
setPlotLayout(length(dataTypes), length(cancerTypes), FALSE)

for (canc in 1:length(cancerTypes)) {
  imputeBoxplot(methyl~impute, canc, ylim=c(0.5,1.5), label="methyl")
  imputeBoxplot(mirna~impute, canc, ylim=c(0.5,1.5), label="mirna")
  imputeBoxplot(mrna~impute, canc, ylim=c(0.5,1.5), label="mrna")
}

labelPlot("Imputation Methods Error", "cancer/imputation method",
          "data type/nrmse")
dev.off()

# 2. Plot each clustering method's accuracy in estimating the clusters
# in the presence of missing data
# -Fifteen box plots (one plot for each cancer/clustering method pair)
# -Three boxes in each plot (one for each imputation method)

clusterBoxplot <- function(formula, canc, clus, ylim, ny=10, label) {
  boxplot(formula=formula, data=cluster,
          subset=cancer==canc&cluster==clus&set==1, ylim=ylim,
          names=imputeTypes)
  grid(nx=NA, ny=ny, col="grey40")
  mtext(paste(cancerTypes[canc], label, sep=", "), side=3, line=-1,
        adj=0.1, cex=0.6, col='grey40')
  
  # Plot survival scores for the intersection of the full data
  if (label %in% c("pval", "ci")) {
    intersectScore <- subset(originalCluster,
                             cancer==canc&cluster==clus&set==2,
                             label)[1,1]
    abline(intersectScore, 0)
    
    unionScores <- subset(originalCluster,
                          cancer==canc&cluster==clus&set==1,
                          label)[, 1]
    lines(unionScores)
    points(unionScores)
  }
}

pdf(file=paste0(plotsFolder, "/cluster.pdf"))
setPlotLayout(4, length(cancerTypes), FALSE)

for (clus in 1:length(clusterTypes)) {
  for (canc in 1:length(cancerTypes)) {
    clusterBoxplot(acc~impute, canc, clus, ylim=c(0,1), label="acc")
    clusterBoxplot(nmi~impute, canc, clus, ylim=c(0,1), label="nmi")
    clusterBoxplot(pval~impute, canc, clus, ylim=c(0,5), label="pval")
    clusterBoxplot(ci~impute, canc, clus, ylim=c(0.45,0.65), label="ci")
  }
  labelPlot(paste("Imputation Methods", clusterTypes[clus], 
                  "Clustering Scores"), "cancer/imputation method",
            "metric/score")
}
dev.off()

# Similarity Plots

similarityBoxplot <- function(formula, canc, ylim, ny=10, label) {
  boxplot(formula=formula, data=similarity,
          subset=cancer==canc&set==1, ylim=ylim,
          names=similarityTypes)
  grid(nx=NA, ny=ny, col="grey40")
  mtext(paste(cancerTypes[canc], label, sep=", "), side=3, line=-1,
        adj=0.1, cex=0.6, col='grey40')
  
  # Plot survival scores for the intersection of the full data
  if (label %in% c("pval", "ci")) {
    intersectScore <- subset(originalSimilarity, cancer==canc&set==2,
                             label)[1,1]
    abline(intersectScore, 0)
    
    unionScores <- subset(originalSimilarity, cancer==canc&set==1,
                          label)[, 1]
    lines(unionScores)
    points(unionScores)
  }
}

pdf(file=paste0(plotsFolder, "/similarity.pdf"))
setPlotLayout(4, length(cancerTypes), FALSE)

for (canc in 1:length(cancerTypes)) {
  similarityBoxplot(acc~similarity, canc, ylim=c(0,1), label="acc")
  similarityBoxplot(nmi~similarity, canc, ylim=c(0,1), label="nmi")
  similarityBoxplot(pval~similarity, canc, ylim=c(0,5), label="pval")
  similarityBoxplot(ci~similarity, canc, ylim=c(0.45,0.65), label="ci")
}

labelPlot("Similarity Methods Clustering Scores",
          "cancer/similarity method", "metric/score")
dev.off()

# 3. Plot the survival curves estimated for each clustering method
# -Six survival plots (two plots for each clustering method)


# Combine the similarity and clustering data sets
# Similarity:
# -add "cluster" column, with "SNF" in each entry
# -rename "similarity" column to "impute"
newSimilarity <- similarity
colnames(newSimilarity)[similarityCols=="similarity"] <- "impute"
newSimilarity$cluster <- "SNF"
newSimilarity$impute <- similarityTypes[newSimilarity$impute]

newCluster <- cluster
newCluster$cluster <- clusterTypes[newCluster$cluster]
newCluster$impute <- imputeTypes[newCluster$impute]

combinedScores <- rbind(newCluster, newSimilarity)
combinedScores <- subset(combinedScores, set==1,
                         c(cancer, impute, seed, cluster, acc, nmi,
                           pval, ci, se.std, pvalcox, coefcox, con_index_ci, con_index_p, bias_corrected_c_index))
combinedScores$method <- interaction(combinedScores$cluster,
                                     combinedScores$impute, drop=TRUE)
methodTypes <- levels(combinedScores$method)
combinedScores$method <- as.numeric(combinedScores$method)

# Remove Inf pval's
combinedScores$pval[is.infinite(combinedScores$pval)] <- NA

write.csv(combinedScores, paste0(results_folder, '/scoresTwoThousandPval.csv'))


# Compare if the mean of method 1 is significantly greater than the
# mean of method 2, for every pair of methods and cancer type
xGreaterY <- function(x, y) {
  t.test(x, y, alternative="greater", paired=TRUE, na.action = na.omit)$p.value < 0.05
}

testScores <- matrix(, 0, 12)

for (canc in 1:length(cancerTypes)) {
  for (met1 in 1:length(methodTypes)) {
    score1 <- subset(combinedScores, cancer==canc&method==met1)
    accVector <- rep.int(0, length(methodTypes))
    nmiVector <- accVector
    pvalVector <- accVector
    ciVector <- accVector
    se.std <- accVector
    pval_cox <- accVector
    coef_cox <- accVector
    con_index_ci <- accVector
    con_index_p <- accVector
    bias_corrected_c_index <- accVector
    for (met2 in (1:length(methodTypes))[-met1]) {
      score2 <- subset(combinedScores, cancer==canc&method==met2)
      accVector[met2] <- xGreaterY(score1$acc, score2$acc)
      nmiVector[met2] <- xGreaterY(score1$nmi, score2$nmi)
      pvalVector[met2] <- xGreaterY(score1$pval, score2$pval)
      ciVector[met2] <- xGreaterY(score1$ci, score2$ci)
      se.std[met2] <- xGreaterY(score1$se.std, score2$se.std)
      pval_cox[met2] <- xGreaterY(score1$pvalcox, score2$pvalcox)
      coef_cox[met2] <- xGreaterY(score1$coefcox, score2$coefcox)
      con_index_ci[met2] <- xGreaterY(score1$con_index_ci, score2$con_index_ci)
      con_index_p[met2] <- xGreaterY(score1$con_index_p, score2$con_index_p)
      bias_corrected_c_index[met2] <- xGreaterY(score1$bias_corrected_c_index, score2$bias_corrected_c_index)
    }
    testScores <- rbind(testScores,
                        c(canc, met1, sum(accVector), sum(nmiVector),
                          sum(pvalVector), sum(ciVector), sum(se.std), sum(pval_cox), sum(coef_cox), sum(con_index_ci),
                          sum(con_index_p), sum(bias_corrected_c_index)))
  }
}

colnames(testScores) <- c("cancer", "method", "acc", "nmi", "pval",
                          "ci")
testScores <- as.data.frame(testScores)
testScores$cancer <- cancerTypes[testScores$cancer]
testScores$method <- methodTypes[testScores$method]

write.table(testScores, "testScores.txt", row.names=FALSE)
