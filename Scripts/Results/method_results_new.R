################################################################################################
# This script will take the results from all folder and make barplots for each strategy and method 
library(ggplot2)
library(reshape2)
library(dplyr)
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
scoresNormalClust <- read.csv(paste0(results_folder, '/scoresTwoThousandClust.csv'))
scoresNormal1000 <- read.csv(paste0(results_folder, '/scoresTwoThousandDup1000.csv'))
scoresNormal3000 <- read.csv(paste0(results_folder, '/scoresTwoThousandDup3000.csv'))


scoresNormalOrig <- read.csv(paste0(results_folder, '/scoresTwoThousandOrig.csv'))
scoresNormalOrigDup <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup.csv'))
scoresNormalOrigDupClust <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDupClust.csv'))

scoresNormalOrigDup1000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup1000.csv'))
scoresNormalOrigDup3000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup3000.csv'))

scoresLUSCOrigDup <- read.csv(paste0(results_folder, '/scoresLUSCOrigDup.csv'))
scoresLUSCNormalDup <- read.csv(paste0(results_folder, '/scoresLUSCNormalDup.csv'))

scoresCombat <- read.csv(paste0(results_folder, '/scoresCombat.csv'))

scoresCombatDup <- read.csv(paste0(results_folder, '/scoresCombatDup.csv'))
scoresCombatOrigDup <- read.csv(paste0(results_folder, '/scoresCombatOrigDup.csv'))

scoresCombatDup3 <- read.csv(paste0(results_folder, '/scoresCombatDup3.csv'))
scoresCombatOrigDup3 <- read.csv(paste0(results_folder, '/scoresCombatOrigDup3.csv'))

scoresCombatDup4 <- read.csv(paste0(results_folder, '/scoresCombatDup4.csv'))
scoresCombatOrigDup4 <- read.csv(paste0(results_folder, '/scoresCombatOrigDup4.csv'))

scoresOrigIntDup <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigIntDup.csv'))

scoresDataTypes <- read.csv(paste0(results_folder, '/scoresDataTypes.csv'))
scoresOrigDataTypes <- read.csv(paste0(results_folder, '/scoresOrigDataTypes.csv'))

scoresCombatDataTypes <- read.csv(paste0(results_folder, '/scoresCombatDataTypes.csv'))
scoresOrigCombatDataTypes <- read.csv(paste0(results_folder, '/scoresOrigCombatDataTypes.csv'))

###################################################################################################
# Get method types
scoresNormal$method <- interaction(scoresNormalClust$cluster,
                                   scoresNormalClust$impute, drop = TRUE)

scoresNormalClust$method <- interaction(scoresNormal$cluster,
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

scoresCombatDup3$method <- interaction(scoresCombatDup3$cluster,
                                      scoresCombatDup3$impute, drop = TRUE)

scoresCombatOrigDup3$method <- interaction(scoresCombatOrigDup3$cluster,
                                          scoresCombatOrigDup3$impute, drop = TRUE)

scoresCombatDup4$method <- interaction(scoresCombatDup4$cluster,
                                       scoresCombatDup4$impute, drop = TRUE)

scoresCombatOrigDup4$method <- interaction(scoresCombatOrigDup4$cluster,
                                           scoresCombatOrigDup4$impute, drop = TRUE)

scoresNormalOrig$method <- interaction(scoresNormalOrig$cluster,
                                       scoresNormalOrig$impute, drop = TRUE)

scoresNormalOrigDup$method <- interaction(scoresNormalOrigDup$cluster,
                                          scoresNormalOrigDup$impute, drop = TRUE)

scoresNormalOrigDupClust$method <- interaction(scoresNormalOrigDupClust$cluster,
                                          scoresNormalOrigDupClust$impute, drop = TRUE)

scoresNormalOrigDup1000$method <- interaction(scoresNormalOrigDup1000$cluster,
                                              scoresNormalOrigDup1000$impute, drop = TRUE)

scoresNormalOrigDup3000$method <- interaction(scoresNormalOrigDup3000$cluster,
                                              scoresNormalOrigDup3000$impute, drop = TRUE)

scoresLUSCOrigDup$method <- interaction(scoresLUSCOrigDup$cluster,
                                        scoresLUSCOrigDup$impute, drop = TRUE)

scoresOrigIntDup$method <- interaction(scoresOrigIntDup$cluster,
                                        scoresOrigIntDup$impute, drop = TRUE)

scoresDataTypes$method <- interaction(scoresDataTypes$cluster,
                                       scoresDataTypes$impute, drop = TRUE)

scoresOrigDataTypes$method <- interaction(scoresOrigDataTypes$cluster,
                                      scoresOrigDataTypes$impute, drop = TRUE)

scoresCombatDataTypes$method <- interaction(scoresCombatDataTypes$cluster,
                                      scoresCombatDataTypes$impute, drop = TRUE)

scoresOrigCombatDataTypes$method <- interaction(scoresOrigCombatDataTypes$cluster,
                                          scoresOrigCombatDataTypes$impute, drop = TRUE)

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

scoresOrigIntDup <- scoresOrigIntDup[complete.cases(scoresOrigIntDup),]

scoresDataTypes <- scoresDataTypes[complete.cases(scoresDataTypes),]
scoresOrigDataTypes <- scoresOrigDataTypes[complete.cases(scoresOrigDataTypes),]

scoresCombatDataTypes <- scoresCombatDataTypes[complete.cases(scoresCombatDataTypes),]
scoresOrigCombatDataTypes <- scoresOrigCombatDataTypes[complete.cases(scoresOrigCombatDataTypes),]


####################################################################################################
# Group by method and get mean for each evaluation method (acc, nmi, pval, ci)

groupbyMethod <- function(data, orig = FALSE, 
                          normal = FALSE,
                          pval = FALSE, 
                          pvalcox = FALSE, 
                          con_index_p = FALSE, 
                          con_index_ci = FALSE, 
                          bias = FALSE, 
                          coef = FALSE,
                          int = FALSE,
                          title) {
  
  if ((orig) && (pval)) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
  
    } else if (orig & int) {
    
      temp <- data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(intAcc, na.rm = T),
                meanNmi = mean(intNmi, na.rm = T))
    
    
      temp_melt <- melt(temp, id.vars = c('method'))
    
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
  
    } else if (orig & pvalcox) {
    
      temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pvalcox, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    
      temp_melt <- melt(temp, id.vars = c('method'))
    
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
  
    } else if (orig & con_index_p) {
   
      temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(con_index_p, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    
      temp_melt <- melt(temp, id.vars = c('method'))
    
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  
    } else if (orig & con_index_ci) {
    
      temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(con_index_ci, na.rm = T))
    
    
      temp_melt <- melt(temp, id.vars = c('method'))
    
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  
    } else if (orig & bias) {
    
      temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(bias_corrected_c_index, na.rm = T))
    
    
      temp_melt <- melt(temp, id.vars = c('method'))
    
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  
    } else if (orig & coef) {
    
      temp <- data %>%
      group_by(method) %>%
      summarise(meanCoef = mean(coefcox, na.rm = T),
                meanStd = mean(se.std, na.rm = T))
    
    
      temp_melt <- melt(temp, id.vars = c('method'))
    
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  
    } else if (normal) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
  }
}


groupbyMethod(scoresNormal, normal = TRUE, title = 'Intersection All Cancers')
groupbyMethod(scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers')
groupbyMethod(scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup')

groupbyMethod(scoresNormal1000, normal = TRUE, title = 'Intersection All Cancers 1000')
groupbyMethod(scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup 1000')

groupbyMethod(scoresNormal3000, normal = TRUE, title = 'Intersection All Cancers 3000')
groupbyMethod(scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup 3000')

groupbyMethod(scoresOrigIntDup, orig = TRUE, int = TRUE, title = 'Union Intersection of Union')

##################################################################################################
# Group by cancer and method and mean of evaluation. 
groupbyCancer <- function(cancer, data, orig = FALSE, 
                          normal = FALSE,
                          pval = FALSE, 
                          pvalcox = FALSE, 
                          con_index_p = FALSE, 
                          con_index_ci = FALSE, 
                          bias = FALSE, 
                          coef = FALSE,
                          int = FALSE,
                          nmi = FALSE, 
                          acc = FALSE,
                          acc_nmi = FALSE,
                          title) {
  
  temp_data <- data[data$cancer == cancer,]
  
  if ((orig) && (pval)) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
    } else if (orig & int) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(intAcc, na.rm = T),
                meanNmi = mean(intNmi, na.rm = T))
    
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
    } else if (orig & acc) {
      
      temp <- temp_data %>%
        group_by(method) %>%
        summarise(meanAcc = mean(intAcc, na.rm = T))
            
      ggplot(data = temp, aes(reorder(method, -meanAcc), meanAcc)) + geom_bar(stat = 'identity') +
        xlab('Method') + ylab('Mean ACC') + ggtitle(title) + theme_538_bar
      
    } else if (orig & nmi) {
      
      temp <- temp_data %>%
        group_by(method) %>%
        summarise(meanNmi = mean(intNmi, na.rm = T))
      
      ggplot(data = temp, aes(reorder(method, -meanNmi), meanNmi)) + geom_bar(stat = 'identity') +
        xlab('Method') + ylab('Mean NMI') + ggtitle(title) + theme_538_bar
    
  } else if (orig & pvalcox) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pvalcox, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  } else if (orig & con_index_p) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanPval = mean(con_index_p, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  } else if (orig & con_index_ci) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(con_index_ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  } else if (orig & bias) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(bias_corrected_c_index, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  } else if (orig & coef) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanCoef = mean(coefcox, na.rm = T),
                meanStd = mean(se.std, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  } else if (normal) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
    
  } else if (normal & pvalcox) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(pvalcox, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
    
    
  } else if (normal & con_index_p) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(con_index_p, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
    
  } else if (normal & con_index_ci) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(pval, na.rm = T),
                meanCi = mean(con_index_ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
  
  } else if (acc_nmi) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T))    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
    
  }
}

groupbyCancer(cancer = 1, scoresNormal, normal = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormal, acc_nmi = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA')
groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA Dup')

groupbyCancer(cancer = 1, scoresNormalClust, normal = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormalClust, acc_nmi = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormalOrigDupClust, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA Dup')


groupbyCancer(cancer = 1, scoresNormal1000, normal = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA Dup 1000')

groupbyCancer(cancer = 1, scoresNormal3000, normal = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA Dup 3000')

groupbyCancer(cancer = 1, scoresOrigIntDup, orig = TRUE, int = TRUE, title = 'BRCA Intersection of Union')
groupbyCancer(cancer = 1, scoresOrigIntDup, orig = TRUE, acc = TRUE, title = 'BRCA Intersection of Union')
groupbyCancer(cancer = 1, scoresOrigIntDup, orig = TRUE, nmi = TRUE, title = 'BRCA Intersection of Union')

## KIRC
groupbyCancer(cancer = 2, scoresNormal, normal = TRUE, title = 'Intersection KIRC')
groupbyCancer(cancer = 2, scoresNormal, acc_nmi = TRUE, title = 'Intersection KIRC')

groupbyCancer(cancer = 2, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union kirc')
groupbyCancer(cancer = 2, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union kirc Dup')

groupbyCancer(cancer = 2, scoresNormalClust, normal = TRUE, title = 'Intersection KIRC')
groupbyCancer(cancer = 2, scoresNormalClust, acc_nmi = TRUE, title = 'Intersection KIRC')

groupbyCancer(cancer = 2, scoresNormalOrigDupClust, orig = TRUE, pval = TRUE, title = 'Pval Union kirc Dup')

groupbyCancer(cancer = 2, scoresNormal1000, normal = TRUE, title = 'Intersection KIRC')
groupbyCancer(cancer = 2, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union kirc Dup 1000')

groupbyCancer(cancer = 2, scoresNormal3000, normal = TRUE, title = 'Intersection KIRC')
groupbyCancer(cancer = 2, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union kirc Dup 3000')

groupbyMethod(scoresCombat, normal = TRUE, title = 'Combat intersection')
groupbyMethod(scoresCombatDup, normal = TRUE, title = 'Combat intersection')
groupbyMethod(scoresCombatOrigDup, orig = TRUE, pval = TRUE, title = 'Combat union duplicate removed')

groupbyMethod(scoresCombatDup3, normal = TRUE, title = 'Combat intersection 3 clusters')
groupbyMethod(scoresCombatOrigDup3, orig = TRUE, pval = TRUE, title = 'Combat union duplicate removed 3 clusters')

groupbyMethod(scoresCombatDup4, normal = TRUE, title = 'Combat intersection 4 clusters')
groupbyMethod(scoresCombatOrigDup4, orig = TRUE, pval = TRUE, title = 'Combat union duplicate removed 4 clusters')

groupbyCancer(cancer = 2, scoresOrigIntDup, orig = TRUE, int = TRUE, title = 'KIRC (Combat) Intersection of Union')
groupbyCancer(cancer = 2, scoresOrigIntDup, orig = TRUE, acc = TRUE, title = 'KIRC (Combat) Intersection of Union')
groupbyCancer(cancer = 2, scoresOrigIntDup, orig = TRUE, nmi = TRUE, title = 'KIRC (Combat) Intersection of Union')

## lihc
groupbyCancer(cancer = 3, scoresNormal, normal = TRUE, title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormal, acc_nmi = TRUE, title = 'Intersection lihc')

groupbyCancer(cancer = 3, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union lihc')
groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union lihc Dup')

groupbyCancer(cancer = 3, scoresNormalClust, normal = TRUE, title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormalClust, acc_nmi = TRUE, title = 'Intersection lihc')

groupbyCancer(cancer = 3, scoresNormalOrigDupClust, orig = TRUE, pval = TRUE, title = 'Pval Union lihc Dup')

groupbyCancer(cancer = 3, scoresNormal1000, normal = TRUE, title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union lihc Dup 1000')

groupbyCancer(cancer = 3, scoresNormal3000, normal = TRUE, title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union lihc Dup 3000')

groupbyCancer(cancer = 3, scoresOrigIntDup, orig = TRUE, int = TRUE, title = 'LIHC Intersection of Union')
groupbyCancer(cancer = 3, scoresOrigIntDup, orig = TRUE, acc = TRUE, title = 'LIHC Intersection of Union')
groupbyCancer(cancer = 3, scoresOrigIntDup, orig = TRUE, nmi = TRUE, title = 'LIHC Intersection of Union')

## luad
groupbyCancer(cancer = 4, scoresNormal, normal = TRUE, title = 'Intersection luad')
groupbyCancer(cancer = 4, scoresNormal, acc_nmi = TRUE, title = 'Intersection luad')

groupbyCancer(cancer = 4, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union luad')
groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union luad Dup')


groupbyCancer(cancer = 4, scoresNormal1000, normal = TRUE, title = 'Intersection luad')
groupbyCancer(cancer = 4, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union luad Dup 1000')

groupbyCancer(cancer = 4, scoresNormal3000, normal = TRUE, title = 'Intersection luad')
groupbyCancer(cancer = 4, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union luad Dup 3000')

groupbyCancer(cancer = 4, scoresOrigIntDup, orig = TRUE, int = TRUE, title = 'LUAD Intersection of Union')
groupbyCancer(cancer = 4, scoresOrigIntDup, orig = TRUE, acc = TRUE, title = 'LUAD Intersection of Union')
groupbyCancer(cancer = 4, scoresOrigIntDup, orig = TRUE, nmi = TRUE, title = 'LUAD Intersection of Union')

# LUSC
groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, pval = TRUE, title = 'Intersection LUSC')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, acc_nmi = TRUE, pval = TRUE, title = 'Intersection LUSC')


groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, pval = TRUE, title = 'Union LUSC duplicates removed')

groupbyCancer(cancer = 5, scoresNormal1000, normal = TRUE, title = 'Intersection lusc')
groupbyCancer(cancer = 5, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union lusc Dup 1000')

groupbyCancer(cancer = 5, scoresNormal3000, normal = TRUE, title = 'Intersection lusc')
groupbyCancer(cancer = 5, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union lusc Dup 3000')

groupbyCancer(cancer = 5, scoresOrigIntDup, orig = TRUE, int = TRUE, title = 'LUSC Intersection of Union')
groupbyCancer(cancer = 5, scoresOrigIntDup, orig = TRUE, acc = TRUE, title = 'LUSC Intersection of Union')
groupbyCancer(cancer = 5, scoresOrigIntDup, orig = TRUE, nmi = TRUE, title = 'LUSC Intersection of Union')

####################################################################################################################
# graphs with error bars

groupError <- function(data, 
                       orig = FALSE, 
                       normal = FALSE,
                       pval = FALSE,
                       all = FALSE,
                       acc_nmi_ci = FALSE,
                       acc_nmi = FALSE,
                       acc = FALSE,
                       nmi = FALSE, 
                       ci = FALSE,
                       title) {
  
  if ((orig) && (pval)) {
    temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  } else if (orig && ci) { 
  
    temp <- data %>%
    group_by(method) %>%
    summarise(meanCi = mean(ci, na.rm = T))
  
    ggplot(data = temp, aes(reorder(method, -meanCi), meanCi)) + geom_bar(stat = 'identity') +
    xlab('Method') + ylab('Mean Ci') + ggtitle(title) + theme_538_bar
  
  } else if (normal && all) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(pval + ci + nmi + acc, na.rm = T),
                stdScore = sd(pval + ci + nmi + acc),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Pval + ci + acc + nmi') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && pval) {
  
  temp <- data %>%
    group_by(method) %>%
    summarise(meanScore = mean(pval , na.rm = T),
              stdScore = sd(pval, na.rm = T),
              count = n())
  
  temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
  temp$low <- temp$meanScore - temp$error
  temp$high <- temp$meanScore + temp$error
  
  
  ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
    xlab('Method') + ylab('Mean Pval') + ggtitle(title) + 
    geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
    theme_538_bar
  
  } else if (normal && acc) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(acc , na.rm = T),
                stdScore = sd(acc, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean ACC') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && ci) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(ci, na.rm = T),
                stdScore = sd(ci, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean CI') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && nmi) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(nmi, na.rm = T),
                stdScore = sd(nmi, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean CI') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && acc_nmi_ci) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(nmi + acc+ ci, na.rm = T),
                stdScore = sd(nmi +acc+ci, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean acc+nmi+ci') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && acc_nmi) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(nmi+acc, na.rm = T),
                stdScore = sd(nmi+acc, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean acc + nmi') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  }
}


groupError(scoresNormal, normal = TRUE, all = TRUE, title = 'Intersection All Cancers')
groupError(scoresNormal, normal = TRUE, acc_nmi_ci = TRUE, title = 'Intersection All Cancers')
groupError(scoresNormal, normal = TRUE, acc_nmi = TRUE, title = 'Intersection All Cancers')
groupError(scoresNormal, normal = TRUE, acc = TRUE, title = 'Intersection All Cancers')
groupError(scoresNormal, normal = TRUE, ci = TRUE, title = 'Intersection All Cancers')
groupError(scoresNormal, normal = TRUE, pval = TRUE, title = 'Intersection All Cancers')
groupError(scoresNormal, normal = TRUE, nmi = TRUE, title = 'Intersection All Cancers')

groupError(scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers')
groupError(scoresNormalOrig, orig = TRUE, ci = TRUE, title = 'Pval Union All Cancers')

groupError(scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup')

groupError(scoresNormal1000, normal = TRUE, title = 'Intersection All Cancers 1000')
groupError(scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup 1000')

groupError(scoresNormal3000, normal = TRUE, title = 'Intersection All Cancers 3000')
groupError(scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup 3000')


##################################################################################################
cancerError <- function(data, 
                        orig = FALSE, 
                        normal = FALSE,
                        pval = FALSE,
                        all = FALSE,
                        acc_nmi_ci = FALSE,
                        acc_nmi = FALSE,
                        acc = FALSE,
                        nmi = FALSE, 
                        ci = FALSE,
                        cancer = NULL,
                        title) {
  
  data <- data[which(data$cancer == cancer),]
  
  if ((orig) && (pval)) {
    temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
    } else if (orig && ci) {
      
      temp <- data %>%
        group_by(method) %>%
        summarise(meanCi = mean(ci, na.rm = T))
      
      ggplot(data = temp, aes(reorder(method, -meanCi), meanCi)) + geom_bar(stat = 'identity') +
        xlab('Method') + ylab('Mean Ci') + ggtitle(title) + theme_538_bar
 
    
    } else if (normal && all) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(pval + ci + nmi + acc, na.rm = T),
                stdScore = sd(pval + ci + nmi + acc),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Pval + ci + acc + nmi') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && pval) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(pval , na.rm = T),
                stdScore = sd(pval, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Pval') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && acc) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(acc , na.rm = T),
                stdScore = sd(acc, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean ACC') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && ci) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(ci, na.rm = T),
                stdScore = sd(ci, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean CI') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && nmi) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(nmi, na.rm = T),
                stdScore = sd(nmi, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean CI') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && acc_nmi_ci) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(nmi + acc+ ci, na.rm = T),
                stdScore = sd(nmi +acc+ci, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean acc+nmi+ci') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  } else if (normal && acc_nmi) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanScore = mean(nmi+acc, na.rm = T),
                stdScore = sd(nmi+acc, na.rm = T),
                count = n())
    
    temp$error <- qt(0.975, df = temp$count - 1)*temp$stdScore/sqrt(temp$count)
    temp$low <- temp$meanScore - temp$error
    temp$high <- temp$meanScore + temp$error
    
    
    ggplot(data = temp, aes(reorder(method, -meanScore), meanScore)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean acc + nmi') + ggtitle(title) + 
      geom_errorbar(mapping = aes(ymin =temp$low, ymax = temp$high), stat = "identity", position = "identity") +
      theme_538_bar
    
  }
}


cancerError(cancer = 1, scoresNormal, normal = TRUE, all = TRUE, title = 'Intersection brca')
cancerError(cancer = 1, scoresNormal, normal = TRUE, acc_nmi_ci = TRUE, title = 'Intersection brca')
cancerError(cancer = 1, scoresNormal, normal = TRUE, acc_nmi = TRUE, title = 'Intersection brca')
cancerError(cancer = 1, scoresNormal, normal = TRUE, acc = TRUE, title = 'Intersection brca')
cancerError(cancer = 1, scoresNormal, normal = TRUE, ci = TRUE, title = 'Intersection brca')
cancerError(cancer = 1, scoresNormal, normal = TRUE, pval = TRUE, title = 'Intersection brca')
cancerError(cancer = 1, scoresNormal, normal = TRUE, nmi = TRUE, title = 'Intersection brca')

cancerError(cancer = 1, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Union brca')
cancerError(cancer = 1, scoresNormalOrig, orig = TRUE, ci = TRUE, title = 'Union brca')
cancerError(cancer = 1, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union brca Dup')
cancerError(cancer = 1, scoresNormalOrigDup, orig = TRUE, ci = TRUE, title = 'Union brca Dup')



cancerError(cancer = 2, scoresNormal, normal = TRUE, all = TRUE, title = 'Intersection kirc')
cancerError(cancer = 2, scoresNormal, normal = TRUE, acc_nmi_ci = TRUE, title = 'Intersection kirc')
cancerError(cancer = 2, scoresNormal, normal = TRUE, acc_nmi = TRUE, title = 'Intersection kirc')
cancerError(cancer = 2, scoresNormal, normal = TRUE, acc = TRUE, title = 'Intersection kirc')
cancerError(cancer = 2, scoresNormal, normal = TRUE, ci = TRUE, title = 'Intersection kirc')
cancerError(cancer = 2, scoresNormal, normal = TRUE, pval = TRUE, title = 'Intersection kirc')
cancerError(cancer = 2, scoresNormal, normal = TRUE, nmi = TRUE, title = 'Intersection kirc')

cancerError(cancer = 2, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Union kirc')
cancerError(cancer = 2, scoresNormalOrig, orig = TRUE, ci = TRUE, title = 'Union kirc')
cancerError(cancer = 2, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union kirc Dup')
cancerError(cancer = 2, scoresNormalOrigDup, orig = TRUE, ci = TRUE, title = 'Union kirc Dup')

groupError(scoresCombat, normal = TRUE, all = TRUE, title = 'intersection combat')
groupError(scoresCombat, normal = TRUE, acc_nmi_ci = TRUE, title = 'intersection combat')
groupError(scoresCombat, normal = TRUE, acc_nmi = TRUE, title = 'intersection combat')
groupError(scoresCombat, normal = TRUE, acc = TRUE, title = 'intersection combat')
groupError(scoresCombat, normal = TRUE, ci = TRUE, title = 'intersection combat')
groupError(scoresCombat, normal = TRUE, pval = TRUE, title = 'intersection combat')
groupError(scoresCombat, normal = TRUE, nmi = TRUE, title = 'intersection combat')

groupError(scoresCombatDup, normal = TRUE, all = TRUE, title = 'intersection combat')
groupError(scoresCombatDup, normal = TRUE, acc_nmi_ci = TRUE, title = 'intersection combat')
groupError(scoresCombatDup, normal = TRUE, acc_nmi = TRUE, title = 'intersection combat')
groupError(scoresCombatDup, normal = TRUE, acc = TRUE, title = 'intersection combat')
groupError(scoresCombatDup, normal = TRUE, ci = TRUE, title = 'intersection combat')
groupError(scoresCombatDup, normal = TRUE, pval = TRUE, title = 'intersection combat')
groupError(scoresCombatDup, normal = TRUE, nmi = TRUE, title = 'intersection combat')

groupError(scoresCombatOrigDup, orig = TRUE, pval = TRUE, title = "Union Combat")
groupError(scoresCombatOrigDup, orig = TRUE, ci = TRUE, title = "Union Combat")


cancerError(cancer = 3, scoresNormal, normal = TRUE, all = TRUE, title = 'Intersection lihc')
cancerError(cancer = 3, scoresNormal, normal = TRUE, acc_nmi_ci = TRUE, title = 'Intersection lihc')
cancerError(cancer = 3, scoresNormal, normal = TRUE, acc_nmi = TRUE, title = 'Intersection lihc')
cancerError(cancer = 3, scoresNormal, normal = TRUE, acc = TRUE, title = 'Intersection lihc')
cancerError(cancer = 3, scoresNormal, normal = TRUE, ci = TRUE, title = 'Intersection lihc')
cancerError(cancer = 3, scoresNormal, normal = TRUE, pval = TRUE, title = 'Intersection lihc')
cancerError(cancer = 3, scoresNormal, normal = TRUE, nmi = TRUE, title = 'Intersection lihc')

cancerError(cancer = 3, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Union lihc')
cancerError(cancer = 3, scoresNormalOrig, orig = TRUE, ci = TRUE, title = 'Union lihc')
cancerError(cancer = 3, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union lihc Dup')
cancerError(cancer = 3, scoresNormalOrigDup, orig = TRUE, ci = TRUE, title = 'Union lihc Dup')



cancerError(cancer = 4, scoresNormal, normal = TRUE, all = TRUE, title = 'Intersection luad')
cancerError(cancer = 4, scoresNormal, normal = TRUE, acc_nmi_ci = TRUE, title = 'Intersection luad')
cancerError(cancer = 4, scoresNormal, normal = TRUE, acc_nmi = TRUE, title = 'Intersection luad')
cancerError(cancer = 4, scoresNormal, normal = TRUE, acc = TRUE, title = 'Intersection luad')
cancerError(cancer = 4, scoresNormal, normal = TRUE, ci = TRUE, title = 'Intersection luad')
cancerError(cancer = 4, scoresNormal, normal = TRUE, pval = TRUE, title = 'Intersection luad')
cancerError(cancer = 4, scoresNormal, normal = TRUE, nmi = TRUE, title = 'Intersection luad')

cancerError(cancer = 4, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Union luad')
cancerError(cancer = 4, scoresNormalOrig, orig = TRUE, ci = TRUE, title = 'Union luad')
cancerError(cancer = 4, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union luad Dup')
cancerError(cancer = 4, scoresNormalOrigDup, orig = TRUE, ci = TRUE, title = 'Union luad Dup')



cancerError(cancer = 5, scoresLUSCNormalDup, normal = TRUE, all = TRUE, title = 'Intersection lusc')
cancerError(cancer = 5, scoresLUSCNormalDup, normal = TRUE, acc_nmi_ci = TRUE, title = 'Intersection lusc')
cancerError(cancer = 5, scoresLUSCNormalDup, normal = TRUE, acc_nmi = TRUE, title = 'Intersection lusc')
cancerError(cancer = 5, scoresLUSCNormalDup, normal = TRUE, acc = TRUE, title = 'Intersection lusc')
cancerError(cancer = 5, scoresLUSCNormalDup, normal = TRUE, ci = TRUE, title = 'Intersection lusc')
cancerError(cancer = 5, scoresLUSCNormalDup, normal = TRUE, pval = TRUE, title = 'Intersection lusc')
cancerError(cancer = 5, scoresLUSCNormalDup, normal = TRUE, nmi = TRUE, title = 'Intersection lusc')

cancerError(cancer = 5, scoresLUSCOrigDup, orig = TRUE, pval = TRUE, title = 'Union luad Dup')
cancerError(cancer = 5, scoresLUSCOrigDup, orig = TRUE, ci = TRUE, title = 'Union luad Dup')

##################################################################################################

# set data types
dataTypes <- c("methyl", "mirna", "mrna")
scoresDataTypes$data_types <- dataTypes[scoresDataTypes$data_types]
scoresOrigDataTypes$data_types <- dataTypes[scoresOrigDataTypes$data_types]
scoresCombatDataTypes$data_types <- dataTypes[scoresCombatDataTypes$data_types]
scoresOrigCombatDataTypes$data_types <- dataTypes[scoresOrigCombatDataTypes$data_types]

# give combat a cancer column 
scoresCombatDataTypes$cancer <- 6
scoresOrigCombatDataTypes$cancer <- 6

# combine scoresDataTypes and combat
scoresAllDataTypes <- rbind(scoresDataTypes, scoresCombatDataTypes)
scoresOrigAllDataTypes <- rbind(scoresOrigDataTypes, scoresOrigCombatDataTypes)


# Group by cancer and method and mean of evaluation. 
groupbyCancer <- function(cancer, 
                          data, 
                          data_types,
                          orig = FALSE, 
                          normal = FALSE,
                          nopval = FALSE, 
                          ci = FALSE, 
                          title) {
  
  temp_data <- data[data$cancer == cancer,]
  temp_data <- temp_data[temp_data$data_types == data_types,]
  
  if (orig && !ci) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
  } else if (orig && ci) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanCi = mean(ci, na.rm = T))
    
    ggplot(data = temp, aes(reorder(method, -meanCi), meanCi)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean CI') + ggtitle(title) + theme_538_bar
    
  } else if (normal && !nopval) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
    
  } else if (normal && nopval) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
    
  }
}

groupbyCancer(cancer = 1, scoresAllDataTypes, data_types = "methyl", normal = TRUE, title = "BRCA Methylation Intersection")
groupbyCancer(cancer = 1, scoresAllDataTypes, data_types = "mirna", normal = TRUE, title = "BRCA mirna Intersecton")
groupbyCancer(cancer = 1, scoresAllDataTypes, data_types = "mrna", normal = TRUE, title = "BRCA mrna Intersection")

groupbyCancer(cancer = 1, scoresAllDataTypes, data_types = "methyl", normal = TRUE, nopval = TRUE, title = "BRCA Methylation Intersection")
groupbyCancer(cancer = 1, scoresAllDataTypes, data_types = "mirna", normal = TRUE, nopval = TRUE, title = "BRCA mirna Intersecton")
groupbyCancer(cancer = 1, scoresAllDataTypes, data_types = "mrna", normal = TRUE, nopval = TRUE, title = "BRCA mrna Intersection")

groupbyCancer(cancer = 1, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, title = "BRCA Methylation Union")
groupbyCancer(cancer = 1, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, title = "BRCA mirna Union")
groupbyCancer(cancer = 1, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, title = "BRCA mrna Union")

groupbyCancer(cancer = 1, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, ci = TRUE, title = "BRCA Methylation Union")
groupbyCancer(cancer = 1, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, ci = TRUE, title = "BRCA mirna Union")
groupbyCancer(cancer = 1, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, ci = TRUE, title = "BRCA mrna Union")

groupbyCancer(cancer = 2, scoresAllDataTypes, data_types = "methyl", normal = TRUE, title = "kirc Methylation Intersection")
groupbyCancer(cancer = 2, scoresAllDataTypes, data_types = "mirna", normal = TRUE, title = "kirc mirna Intersecton")
groupbyCancer(cancer = 2, scoresAllDataTypes, data_types = "mrna", normal = TRUE, title = "kirc mrna Intersection")

groupbyCancer(cancer = 2, scoresAllDataTypes, data_types = "methyl", normal = TRUE, nopval = TRUE, title = "kirc Methylation Intersection")
groupbyCancer(cancer = 2, scoresAllDataTypes, data_types = "mirna", normal = TRUE, nopval = TRUE, title = "kirc mirna Intersecton")
groupbyCancer(cancer = 2, scoresAllDataTypes, data_types = "mrna", normal = TRUE, nopval = TRUE, title = "kirc mrna Intersection")

groupbyCancer(cancer = 2, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, title = "kirc Methylation Union")
groupbyCancer(cancer = 2, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, title = "kirc mirna Union")
groupbyCancer(cancer = 2, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, title = "kirc mrna Union")

groupbyCancer(cancer = 2, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, ci = TRUE, title = "kirc Methylation Union")
groupbyCancer(cancer = 2, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, ci = TRUE, title = "kirc mirna Union")
groupbyCancer(cancer = 2, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, ci = TRUE, title = "kirc mrna Union")

groupbyCancer(cancer = 6, scoresAllDataTypes, data_types = "methyl", normal = TRUE, title = "combat Methylation Intersection")
groupbyCancer(cancer = 6, scoresAllDataTypes, data_types = "mirna", normal = TRUE, title = "combat mirna Intersecton")
groupbyCancer(cancer = 6, scoresAllDataTypes, data_types = "mrna", normal = TRUE, title = "combat mrna Intersection")

groupbyCancer(cancer = 6, scoresAllDataTypes, data_types = "methyl", normal = TRUE, nopval = TRUE, title = "combat Methylation Intersection")
groupbyCancer(cancer = 6, scoresAllDataTypes, data_types = "mirna", normal = TRUE, nopval = TRUE, title = "combat mirna Intersecton")
groupbyCancer(cancer = 6, scoresAllDataTypes, data_types = "mrna", normal = TRUE, nopval = TRUE, title = "combat mrna Intersection")

groupbyCancer(cancer = 6, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, title = "combat Methylation Union")
groupbyCancer(cancer = 6, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, title = "combat mirna Union")
groupbyCancer(cancer = 6, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, title = "combat mrna Union")

groupbyCancer(cancer = 6, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, ci = TRUE, title = "combat Methylation Union")
groupbyCancer(cancer = 6, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, ci = TRUE, title = "combat mirna Union")
groupbyCancer(cancer = 6, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, ci = TRUE, title = "combat mrna Union")

groupbyCancer(cancer = 3, scoresAllDataTypes, data_types = "methyl", normal = TRUE, title = "lihc Methylation Intersection")
groupbyCancer(cancer = 3, scoresAllDataTypes, data_types = "mirna", normal = TRUE, title = "lihc mirna Intersecton")
groupbyCancer(cancer = 3, scoresAllDataTypes, data_types = "mrna", normal = TRUE, title = "lihc mrna Intersection")

groupbyCancer(cancer = 3, scoresAllDataTypes, data_types = "methyl", normal = TRUE, nopval = TRUE, title = "lihc Methylation Intersection")
groupbyCancer(cancer = 3, scoresAllDataTypes, data_types = "mirna", normal = TRUE, nopval = TRUE, title = "lihc mirna Intersecton")
groupbyCancer(cancer = 3, scoresAllDataTypes, data_types = "mrna", normal = TRUE, nopval = TRUE, title = "lihc mrna Intersection")

groupbyCancer(cancer = 3, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, title = "lihc Methylation Union")
groupbyCancer(cancer = 3, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, title = "lihc mirna Union")
groupbyCancer(cancer = 3, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, title = "lihc mrna Union")

groupbyCancer(cancer = 3, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, ci = TRUE, title = "lihc Methylation Union")
groupbyCancer(cancer = 3, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, ci = TRUE, title = "lihc mirna Union")
groupbyCancer(cancer = 3, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, ci = TRUE, title = "lihc mrna Union")

groupbyCancer(cancer = 4, scoresAllDataTypes, data_types = "methyl", normal = TRUE, title = "luad Methylation Intersection")
groupbyCancer(cancer = 4, scoresAllDataTypes, data_types = "mirna", normal = TRUE, title = "luad mirna Intersecton")
groupbyCancer(cancer = 4, scoresAllDataTypes, data_types = "mrna", normal = TRUE, title = "luad mrna Intersection")

groupbyCancer(cancer = 4, scoresAllDataTypes, data_types = "methyl", normal = TRUE, nopval = TRUE,title = "luad Methylation Intersection")
groupbyCancer(cancer = 4, scoresAllDataTypes, data_types = "mirna", normal = TRUE, nopval = TRUE, title = "luad mirna Intersecton")
groupbyCancer(cancer = 4, scoresAllDataTypes, data_types = "mrna", normal = TRUE, nopval = TRUE, title = "luad mrna Intersection")

groupbyCancer(cancer = 4, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, title = "luad Methylation Union")
groupbyCancer(cancer = 4, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, title = "luad mirna Union")
groupbyCancer(cancer = 4, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, title = "luad mrna Union")

groupbyCancer(cancer = 4, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, ci = TRUE, title = "luad Methylation Union")
groupbyCancer(cancer = 4, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, ci = TRUE,  title = "luad mirna Union")
groupbyCancer(cancer = 4, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, ci = TRUE, title = "luad mrna Union")

groupbyCancer(cancer = 5, scoresAllDataTypes, data_types = "methyl", normal = TRUE, title = "lusc Methylation Intersection")
groupbyCancer(cancer = 5, scoresAllDataTypes, data_types = "mirna", normal = TRUE, title = "lusc mirna Intersecton")
groupbyCancer(cancer = 5, scoresAllDataTypes, data_types = "mrna", normal = TRUE, title = "lusc mrna Intersection")

groupbyCancer(cancer = 5, scoresAllDataTypes, data_types = "methyl", normal = TRUE, nopval = TRUE, title = "lusc Methylation Intersection")
groupbyCancer(cancer = 5, scoresAllDataTypes, data_types = "mirna", normal = TRUE, nopval = TRUE, title = "lusc mirna Intersecton")
groupbyCancer(cancer = 5, scoresAllDataTypes, data_types = "mrna", normal = TRUE, nopval = TRUE, title = "lusc mrna Intersection")

groupbyCancer(cancer = 5, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, title = "lusc Methylation Union")
groupbyCancer(cancer = 5, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, title = "lusc mirna Union")
groupbyCancer(cancer = 5, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, title = "lusc mrna Union")

groupbyCancer(cancer = 5, scoresOrigAllDataTypes, data_types = "methyl", orig = TRUE, ci = TRUE, title = "lusc Methylation Union")
groupbyCancer(cancer = 5, scoresOrigAllDataTypes, data_types = "mirna", orig = TRUE, ci = TRUE, title = "lusc mirna Union")
groupbyCancer(cancer = 5, scoresOrigAllDataTypes, data_types = "mrna", orig = TRUE, ci = TRUE, title = "lusc mrna Union")

