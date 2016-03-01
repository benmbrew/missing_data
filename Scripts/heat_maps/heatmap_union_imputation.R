######### PCA of intersected data with randomly removed and imputed.# Initialize folders
library(impute)
library(gplots)
# Initialize folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features",
                    "evaluate_original_imputation", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
resultsFolder <- paste(testFolder, "Results", sep="/")
setwd(dataFolder)
load('imputedDataPCA.RData')


# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
data_types <- c("methyl", "mirna", "mrna")
num_views <- 3


source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))


heatMapIndividual <- function(data, name, name2){ 
  data <- t(data)
  data_heatmap <- heatmap.2(data, 
                            dendrogram="none", 
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

heatMapIndividual(brca_knn_methyl, name = 'BRCA KNN Methyl', name2 = 'methylation')
heatMapIndividual(brca_knn_mirna, name = 'BRCA KNN mirna', name2 = 'mirna')
heatMapIndividual(brca_knn_mrna, name = 'BRCA KNN mrna', name2 = 'mrna')

heatMapIndividual(brca_lls_methyl, name = 'BRCA lls Methyl', name2 = 'methylation')
heatMapIndividual(brca_lls_mirna, name = 'BRCA lls mirna', name2 = 'mirna')
heatMapIndividual(brca_lls_mrna, name = 'BRCA lls mrna', name2 = 'mrna')

heatMapIndividual(brca_rand_methyl, name = 'BRCA rand Methyl', name2 = 'methylation')
heatMapIndividual(brca_rand_mirna, name = 'BRCA rand mirna', name2 = 'mirna')
heatMapIndividual(brca_rand_mrna, name = 'BRCA rand mrna', name2 = 'mrna')

heatMapIndividual(kirc_knn_methyl, name = 'kirc KNN Methyl', name2 = 'methylation')
heatMapIndividual(kirc_knn_mirna, name = 'kirc KNN mirna', name2 = 'mirna')
heatMapIndividual(kirc_knn_mrna, name = 'kirc KNN mrna', name2 = 'mrna')

heatMapIndividual(kirc_lls_methyl, name = 'kirc lls Methyl', name2 = 'methylation')
heatMapIndividual(kirc_lls_mirna, name = 'kirc lls mirna', name2 = 'mirna')
heatMapIndividual(kirc_lls_mrna, name = 'kirc lls mrna', name2 = 'mrna')

heatMapIndividual(kirc_rand_methyl, name = 'kirc rand Methyl', name2 = 'methylation')
heatMapIndividual(kirc_rand_mirna, name = 'kirc rand mirna', name2 = 'mirna')
heatMapIndividual(kirc_rand_mrna, name = 'kirc rand mrna', name2 = 'mrna')

heatMapIndividual(lihc_knn_methyl, name = 'lihc KNN Methyl', name2 = 'methylation')
heatMapIndividual(lihc_knn_mirna, name = 'lihc KNN mirna', name2 = 'mirna')
heatMapIndividual(lihc_knn_mrna, name = 'lihc KNN mrna', name2 = 'mrna')

heatMapIndividual(lihc_lls_methyl, name = 'lihc lls Methyl', name2 = 'methylation')
heatMapIndividual(lihc_lls_mirna, name = 'lihc lls mirna', name2 = 'mirna')
heatMapIndividual(lihc_lls_mrna, name = 'lihc lls mrna', name2 = 'mrna')

heatMapIndividual(lihc_rand_methyl, name = 'lihc rand Methyl', name2 = 'methylation')
heatMapIndividual(lihc_rand_mirna, name = 'lihc rand mirna', name2 = 'mirna')
heatMapIndividual(lihc_rand_mrna, name = 'lihc rand mrna', name2 = 'mrna')

heatMapIndividual(luad_knn_methyl, name = 'luad KNN Methyl', name2 = 'methylation')
heatMapIndividual(luad_knn_mirna, name = 'luad KNN mirna', name2 = 'mirna')
heatMapIndividual(luad_knn_mrna, name = 'luad KNN mrna', name2 = 'mrna')

heatMapIndividual(luad_lls_methyl, name = 'luad lls Methyl', name2 = 'methylation')
heatMapIndividual(luad_lls_mirna, name = 'luad lls mirna', name2 = 'mirna')
heatMapIndividual(luad_lls_mrna, name = 'luad lls mrna', name2 = 'mrna')

heatMapIndividual(luad_rand_methyl, name = 'luad rand Methyl', name2 = 'methylation')
heatMapIndividual(luad_rand_mirna, name = 'luad rand mirna', name2 = 'mirna')
heatMapIndividual(luad_rand_mrna, name = 'luad rand mrna', name2 = 'mrna')

