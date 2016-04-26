################################################################################################
# This script will take the results from all folder and make barplots for each strategy and method 
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
################################################################################################
# Initialize folders, 
home_folder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
project_folder <- paste(home_folder, 'Projects/SNF/NM_2015', sep = '/')
scripts_folder <- paste(project_folder, "Scripts", '06_Two_Thousand_Features', sep = "/")
plots_folder <- paste(scripts_folder, "analyze_results/Plots", sep ="/")
results_folder <- paste(project_folder, 'Scripts/06_Results', sep = '/')
source(paste0(results_folder, '/Lib/helpers.R'))


# Load data
scoresNormal <- read.csv(paste0(results_folder, '/scoresTwoThousand.csv'))
scoresNormal1000 <- read.csv(paste0(results_folder, '/scoresTwoThousandDup1000.csv'))
scoresNormal3000 <- read.csv(paste0(results_folder, '/scoresTwoThousandDup3000.csv'))


scoresNormalOrig <- read.csv(paste0(results_folder, '/scoresTwoThousandOrig.csv'))
scoresNormalOrigDup <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup.csv'))
scoresNormalOrigDup1000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup1000.csv'))
scoresNormalOrigDup3000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup3000.csv'))

scoresLUSCOrigDup <- read.csv(paste0(results_folder, '/scoresLUSCOrigDup.csv'))
scoresLUSCNormalDup <- read.csv(paste0(results_folder, '/scoresLUSCNormalDup.csv'))

scoresCombat <- read.csv(paste0(results_folder, '/scoresCombat.csv'))

scoresCombatDup <- read.csv(paste0(results_folder, '/scoresCombatDup.csv'))
scoresCombatOrigDup <- read.csv(paste0(results_folder, '/scoresCombatOrigDup.csv'))

###################################################################################################
# Get method types
scoresNormal$method <- interaction(scoresNormal$cluster,
                                   scoresNormal$impute, drop = TRUE)

scoresNormal1000$method <- interaction(scoresNormal1000$cluster,
                                       scoresNormal1000$impute, drop = TRUE)

scoresNormal3000$method <- interaction(scoresNormal3000$cluster,
                                       scoresNormal3000$impute, drop = TRUE)

scoresLUSCNormalDup$method <- interaction(scoresLUSCNormalDup$cluster,
                                          scoresLUSCNormalDup$impute, drop = TRUE)

scoresCombat$method <- interaction(scoresCombat$cluster,
                                   scoresCombat$impute, drop = TRUE)

scoresCombatDup$method <- interaction(scoresCombatDup$cluster,
                                      scoresCombatDup$impute, drop = TRUE)

scoresCombatOrigDup$method <- interaction(scoresCombatOrigDup$cluster,
                                          scoresCombatOrigDup$impute, drop = TRUE)


scoresNormalOrig$method <- interaction(scoresNormalOrig$cluster,
                                       scoresNormalOrig$impute, drop = TRUE)

scoresNormalOrigDup$method <- interaction(scoresNormalOrigDup$cluster,
                                          scoresNormalOrigDup$impute, drop = TRUE)

scoresNormalOrigDup1000$method <- interaction(scoresNormalOrigDup1000$cluster,
                                              scoresNormalOrigDup1000$impute, drop = TRUE)

scoresNormalOrigDup3000$method <- interaction(scoresNormalOrigDup3000$cluster,
                                              scoresNormalOrigDup3000$impute, drop = TRUE)

scoresLUSCOrigDup$method <- interaction(scoresLUSCOrigDup$cluster,
                                        scoresLUSCOrigDup$impute, drop = TRUE)

####################################################################################################
# remove NAs 
scoresNormal <- scoresNormal[complete.cases(scoresNormal),]
scoresNormal1000 <- scoresNormal1000[complete.cases(scoresNormal1000),]
scoresNormal3000 <- scoresNormal3000[complete.cases(scoresNormal3000),]

scoresNormalOrig <- scoresNormalOrig[complete.cases(scoresNormalOrig),]
scoresNormalOrigDup <- scoresNormalOrigDup[complete.cases(scoresNormalOrigDup),]
scoresNormalOrigDup1000 <- scoresNormalOrigDup1000[complete.cases(scoresNormalOrigDup1000),]
scoresNormalOrigDup3000 <- scoresNormalOrigDup3000[complete.cases(scoresNormalOrigDup3000),]

scoresLUSCNormalDup <- scoresLUSCNormalDup[complete.cases(scoresLUSCNormalDup),]
scoresLUSCOrigDup <- scoresLUSCOrigDup[complete.cases(scoresLUSCOrigDup),]

scoresCombat <- scoresCombat[complete.cases(scoresCombat),]
scoresCombatDup <- scoresCombatDup[complete.cases(scoresCombatDup),]
scoresCombatOrigDup <- scoresCombatOrigDup[complete.cases(scoresCombatOrigDup),]

####################################################################################################
# For each cancer, compare all scores together 

cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")

isGreater <- function(data, cancer) {
  
  cancer <- data[which(data$cancer == cancer),]
  
  temp <- gather()
  
  temp <- cancer %>%
    group_by(method) %>%
    summarise(meanScore = mean(acc + nmi + pval + ci, na.rm = T))
  
  
}
for (i in 1:length(cancerTypes)) {
  
}