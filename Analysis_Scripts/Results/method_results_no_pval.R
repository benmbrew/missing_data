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
scoresNormalOrig <- read.csv(paste0(results_folder, '/scoresTwoThousandOrig.csv'))
scoresNormalOrigInt <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigInt.csv'))
scoresNormalOrigIntDup <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigIntDup.csv'))
scoresNormalOrigClust <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigClust.csv'))
scoresNormalOrigDup <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup.csv'))
scoresNormalOrigDup1000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup1000.csv'))
scoresNormalOrigDup3000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup3000.csv'))

scoresNormalOrigNA <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandNA.csv'))
scoresLUSCOrigDup <- read.csv(paste0(results_folder, '/scoresLUSCOrigDup.csv'))
scoresLUSCNormalDup <- read.csv(paste0(results_folder, '/scoresLUSCNormalDup.csv'))
scoresNormalOrigDupSeed <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandSeed.csv'))



scoresCombat <- read.csv(paste0(results_folder, '/scoresCombat.csv'))
scoresCombatOrig <- read.csv(paste0(results_folder, '/scoresCombatOrig.csv'))
scoresCombatOrigDup <- read.csv(paste0(results_folder, '/scoresCombatOrigDup.csv'))
scoresCombatOrigDupSeed <- read.csv(paste0(results_folder, '/scoresCombatOrigDupSeed.csv'))



scoresGender <- read.csv(paste0(results_folder, '/scoresGender.csv'))
scoresGenderOrig <- read.csv(paste0(results_folder, '/scoresGenderOrig.csv'))

# similarity with 300 features
scoresUnion3000 <- read.csv(paste0(results_folder, '/scoresunion300.csv'))
scoresNormal3000 <- read.csv(paste0(results_folder, '/scoresCom300.csv'))


###################################################################################################
# Get method types
scoresNormal$method <- interaction(scoresNormal$cluster,
                                   scoresNormal$impute, drop = TRUE)

scoresNormalOrig$method <- interaction(scoresNormalOrig$cluster,
                                   scoresNormalOrig$impute, drop = TRUE)

scoresNormalOrigInt$method <- interaction(scoresNormalOrigInt$cluster,
                                       scoresNormalOrigInt$impute, drop = TRUE)

scoresNormalOrigIntDup$method <- interaction(scoresNormalOrigIntDup$cluster,
                                          scoresNormalOrigIntDup$impute, drop = TRUE)


scoresNormalOrigClust$method <- interaction(scoresNormalOrigClust$cluster,
                                       scoresNormalOrigClust$impute, drop = TRUE)

scoresNormalOrigDup$method <- interaction(scoresNormalOrigDup$cluster,
                                            scoresNormalOrigDup$impute, drop = TRUE)

scoresNormalOrigDup1000$method <- interaction(scoresNormalOrigDup1000$cluster,
                                          scoresNormalOrigDup1000$impute, drop = TRUE)

scoresNormalOrigDup3000$method <- interaction(scoresNormalOrigDup3000$cluster,
                                          scoresNormalOrigDup3000$impute, drop = TRUE)


scoresNormalOrigNA$method <- interaction(scoresNormalOrigNA$cluster,
                                          scoresNormalOrigNA$impute, drop = TRUE)


scoresLUSCOrigDup$method <- interaction(scoresLUSCOrigDup$cluster,
                                 scoresLUSCOrigDup$impute, drop = TRUE)

scoresLUSCNormalDup$method <- interaction(scoresLUSCNormalDup$cluster,
                                        scoresLUSCNormalDup$impute, drop = TRUE)


scoresNormalOrigDupSeed$method <- interaction(scoresNormalOrigDupSeed$cluster,
                                       scoresNormalOrigDupSeed$impute, drop = TRUE)

scoresCombat$method <- interaction(scoresCombat$cluster,
                            scoresCombat$impute, drop = TRUE)

scoresCombatOrig$method <- interaction(scoresCombatOrig$cluster,
                                scoresCombatOrig$impute, drop = TRUE)

scoresCombatOrigDup$method <- interaction(scoresCombatOrigDup$cluster,
                                   scoresCombatOrigDup$impute, drop = TRUE)

scoresCombatOrigDupSeed$method <- interaction(scoresCombatOrigDupSeed$cluster,
                                       scoresCombatOrigDupSeed$impute, drop = TRUE)

scoresGender$method <- interaction(scoresGender$cluster,
                                   scoresGender$impute, drop = TRUE)

scoresGenderOrig$method <- interaction(scoresGenderOrig$cluster,
                                       scoresGenderOrig$impute, drop = TRUE)

# separate into male and female for complete data
scoresGenderMale <- scoresGender[scoresGender$gender == 1,]
scoresGenderFemale <- scoresGender[scoresGender$gender == 2,]

# separate into male and female
scoresGenderOrigMale <- scoresGenderOrig[scoresGenderOrig$gender == 1,]
scoresGenderOrigFemale <- scoresGenderOrig[scoresGenderOrig$gender == 2,]

scoresNormal3000$method <- interaction(scoresNormal3000$cluster,
                                      scoresNormal3000$impute, drop = TRUE)

scoresUnion3000$method <- interaction(scoresUnion3000$cluster,
                                       scoresUnion3000$impute, drop = TRUE)


###################################################################################################
# Drop acc and nmi from original data, as it's not needed for a evaluation
scoresNormalOrig$acc <- NULL
scoresNormalOrig$nmi <- NULL
scoresNormalOrigClust$acc <- NULL
scoresNormalOrigClust$nmi <- NULL
scoresCombatOrig$acc <- NULL
scoresCombatOrig$nmi <- NULL
scoresGenderOrig$acc <- NULL
scoresGenderOrig$nmi <- NULL
scoresNormalOrigDup$acc <- NULL
scoresNormalOrigDup1000$nmi <- NULL
scoresNormalOrigDup1000$acc <- NULL
scoresNormalOrigDup3000$nmi <- NULL
scoresNormalOrigDup3000$acc <- NULL
scoresNormalOrigDup$nmi <- NULL
scoresNormalOrigNA$acc <- NULL
scoresNormalOrigNA$nmi <- NULL
scoresUnion3000$acc <- NULL
scoresUnion3000$nmi <- NULL

####################################################################################################
# Group by method and get mean for each evaluation method (acc, nmi, pval, ci)

groupbyMethod <- function(data, orig = FALSE, 
                          pval = FALSE,
                          normal = FALSE,
                          acc_nmi = FALSE,
                          acc = FALSE,
                          nmi = FALSE,
                          ci = FALSE,
                          title) {
  
   if ((orig) && (pval)) {
      temp <- data %>%
        group_by(method) %>%
        summarise(meanCi = mean(ci, na.rm = T))
      
      temp_melt <- melt(temp, id.vars = c('method'))
      
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
         xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
     
             
  } else  if (normal) {
              
    temp <- data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
                
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    xlab('Method') + ggtitle(title) + theme_538_bar
  
    } else if (acc_nmi) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T))    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
  
   } else if (acc) {
  
  temp <- data %>%
    group_by(method) %>%
    summarise(meanAcc = mean(acc, na.rm = T))
    
  temp_melt <- melt(temp, id.vars = c('method'))
  
  ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    xlab('Method') + ggtitle(title) + theme_538_bar


  } else if (nmi) {
  
  temp <- data %>%
    group_by(method) %>%
    summarise(meanNmi = mean(nmi, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('method'))
  
  ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    xlab('Method') + ggtitle(title) + theme_538_bar
  
  
  } else if (ci) {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
  }


}


groupbyMethod(scoresNormal, normal = TRUE, title = 'Intersection All Cancers')
groupbyMethod(scoresNormal, acc_nmi = TRUE, title = 'Intersection All Cancers')
groupbyMethod(scoresNormal, acc = TRUE, title = 'Intersection All Cancers')
groupbyMethod(scoresNormal, nmi = TRUE, title = 'Intersection All Cancers')


groupbyMethod(scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup')

##################################################################################################
# Group by cancer and method and mean of evaluation. 
groupbyCancer <- function(cancer, 
                          data, orig = FALSE, 
                          pval = FALSE,
                          normal = FALSE,
                          acc_nmi = FALSE,
                          acc = FALSE,
                          nmi = FALSE,
                          ci = FALSE,
                          title) {
  
  temp_data <- data[data$cancer == cancer,]
  
  if ((orig) && (pval)) {
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
    
    
  } else  if (normal) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
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
    
  } else if (acc) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
    
    
  } else if (nmi) {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanNmi = mean(nmi, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar  
  
  } else if (ci) {
  
  temp <- temp_data %>%
    group_by(method) %>%
    summarise(meanCi = mean(ci, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('method'))
  
  ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    xlab('Method') + ggtitle(title) + theme_538_bar
  }
  
  
}

groupbyCancer(cancer = 1, scoresNormal, normal = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormal, acc_nmi = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormal, acc = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormal, nmi = TRUE, title = 'Intersection BRCA')
groupbyCancer(cancer = 1, scoresNormal, ci = TRUE, title = 'Intersection BRCA')


groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union BRCA Dup')

## KIRC
groupbyCancer(cancer = 2, scoresNormal, normal = TRUE, pval = TRUE, title = 'Intersection kirc')
groupbyCancer(cancer = 2, scoresNormal, acc_nmi = TRUE, pval = TRUE, title = 'Intersection kirc')
groupbyCancer(cancer = 2, scoresNormal, acc = TRUE, pval = TRUE, title = 'Intersection kirc')
groupbyCancer(cancer = 2, scoresNormal, nmi = TRUE, pval = TRUE, title = 'Intersection kirc')
groupbyCancer(cancer = 2, scoresNormal, ci = TRUE, pval = TRUE, title = 'Intersection kirc')

groupbyCancer(cancer = 2, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union kirc Dup')

# combat
groupbyMethod(scoresCombat, normal = TRUE, title = 'Combat intersection')
groupbyMethod(scoresCombat, acc_nmi = TRUE, title = 'Combat intersection')
groupbyMethod(scoresCombat, acc = TRUE, title = 'Combat intersection')
groupbyMethod(scoresCombat, nmi = TRUE, title = 'Combat intersection')
groupbyMethod(scoresCombat, ci = TRUE,  title = 'Combat intersection')

groupbyMethod(scoresCombatOrigDup, orig = TRUE, pval = TRUE, title = 'Combat union duplicate removed')

## lihc
groupbyCancer(cancer = 3, scoresNormal, normal = TRUE,  title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormal, acc_nmi = TRUE,  title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormal, acc = TRUE,  title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormal, nmi = TRUE,  title = 'Intersection lihc')
groupbyCancer(cancer = 3, scoresNormal, ci = TRUE,  title = 'Intersection lihc')


groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union lihc Dup')

## luad
groupbyCancer(cancer = 4, scoresNormal, normal = TRUE, title = 'Intersection luad')
groupbyCancer(cancer = 4, scoresNormal, acc_nmi = TRUE, title = 'Intersection luad')
groupbyCancer(cancer = 4, scoresNormal, acc = TRUE, title = 'Intersection luad')
groupbyCancer(cancer = 4, scoresNormal, nmi = TRUE, title = 'Intersection luad')
groupbyCancer(cancer = 4, scoresNormal, ci = TRUE, title = 'Intersection luad')


groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Union luad Dup')

# LUSC
groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, title = 'Intersection LUSC')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, acc_nmi = TRUE, title = 'Intersection LUSC')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, acc = TRUE, title = 'Intersection LUSC')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, nmi = TRUE, title = 'Intersection LUSC')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, ci = TRUE, title = 'Intersection LUSC')

groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, pval = TRUE, title = 'Union LUSC')

###############################################################################################################

## By Cancer
groupByCancerInt <- function(cancer, title) {
  temp.data <- scoresNormalOrigIntDup[scoresNormalOrigIntDup$cancer == cancer,]
  
  temp <- temp.data %>%
    group_by(method) %>%
    summarise(meanAcc = mean(intAcc, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('method'))
  
  ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title) + theme_538_bar
  
}

groupByCancerInt(1, title = 'BRCA Intersection of Union')
groupByCancerInt(2, title = 'KIRC Intersection of Union')
groupByCancerInt(3, title = 'LIHC Intersection of Union')
groupByCancerInt(4, title = 'LUAD Intersection of Union')
groupByCancerInt(5, title = 'LUSC Intersection of Union')


## By Cancer
groupByCancerInt <- function(cancer, title) {
  temp.data <- scoresNormalOrigIntDup[scoresNormalOrigIntDup$cancer == cancer,]
  
  temp <- temp.data %>%
    group_by(method) %>%
    summarise(meanNmi = mean(intNmi, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('method'))
  
  ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title) + theme_538_bar
  
}

groupByCancerInt(1, title = 'BRCA Intersection of Union')
groupByCancerInt(2, title = 'KIRC Intersection of Union')
groupByCancerInt(3, title = 'LIHC Intersection of Union')
groupByCancerInt(4, title = 'LUAD Intersection of Union')
groupByCancerInt(5, title = 'LUSC Intersection of Union')



