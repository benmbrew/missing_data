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
scoresNormal <- read.csv(paste0(results_folder, '/scoresTwoThousandDup1000.csv'))
scoresNormal <- read.csv(paste0(results_folder, '/scoresTwoThousandDup3000.csv'))


scoresNormalOrig <- read.csv(paste0(results_folder, '/scoresTwoThousandOrig.csv'))
# scoresNormalOrigInt <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigInt.csv'))
scoresNormalOrigIntDup <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigIntDup.csv'))
# scoresNormalOrigClust <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigClust.csv'))
scoresNormalOrigDup <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup.csv'))
scoresNormalOrigDup1000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup1000.csv'))
scoresNormalOrigDup3000 <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDup3000.csv'))

# scoresNormalOrigNA <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandNA.csv'))
scoresLUSCOrigDup <- read.csv(paste0(results_folder, '/scoresLUSCOrigDup.csv'))
scoresLUSCNormalDup <- read.csv(paste0(results_folder, '/scoresLUSCNormalDup.csv'))
# scoresNormalOrigDupSeed <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandSeed.csv'))



scoresCombat <- read.csv(paste0(results_folder, '/scoresCombat.csv'))
scoresCombatDup <- read.csv(paste0(results_folder, '/scoresCombatDup.csv'))
scoresCombatOrig <- read.csv(paste0(results_folder, '/scoresCombatOrig.csv'))
# scoresCombatOrigDupSeed <- read.csv(paste0(results_folder, '/scoresCombatOrigDupSeed.csv'))



# scoresGender <- read.csv(paste0(results_folder, '/scoresGender.csv'))
# scoresGenderOrig <- read.csv(paste0(results_folder, '/scoresGenderOrig.csv'))

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


# scoresNormalOrigClust$method <- interaction(scoresNormalOrigClust$cluster,
#                                        scoresNormalOrigClust$impute, drop = TRUE)

scoresNormalOrigDup$method <- interaction(scoresNormalOrigDup$cluster,
                                            scoresNormalOrigDup$impute, drop = TRUE)

scoresNormalOrigDup1000$method <- interaction(scoresNormalOrigDup1000$cluster,
                                          scoresNormalOrigDup1000$impute, drop = TRUE)

scoresNormalOrigDup3000$method <- interaction(scoresNormalOrigDup3000$cluster,
                                          scoresNormalOrigDup3000$impute, drop = TRUE)


# scoresNormalOrigNA$method <- interaction(scoresNormalOrigNA$cluster,
#                                           scoresNormalOrigNA$impute, drop = TRUE)


scoresLUSCOrigDup$method <- interaction(scoresLUSCOrigDup$cluster,
                                 scoresLUSCOrigDup$impute, drop = TRUE)

scoresLUSCNormalDup$method <- interaction(scoresLUSCNormalDup$cluster,
                                        scoresLUSCNormalDup$impute, drop = TRUE)


# scoresNormalOrigDupSeed$method <- interaction(scoresNormalOrigDupSeed$cluster,
#                                        scoresNormalOrigDupSeed$impute, drop = TRUE)

scoresCombat$method <- interaction(scoresCombat$cluster,
                            scoresCombat$impute, drop = TRUE)

scoresCombatDup$method <- interaction(scoresCombatDup$cluster,
                                   scoresCombatDup$impute, drop = TRUE)

# scoresCombatOrig$method <- interaction(scoresCombatOrig$cluster,
#                                 scoresCombatOrig$impute, drop = TRUE)

scoresCombatOrigDup$method <- interaction(scoresCombatOrigDup$cluster,
                                   scoresCombatOrigDup$impute, drop = TRUE)

# scoresCombatOrigDupSeed$method <- interaction(scoresCombatOrigDupSeed$cluster,
#                                        scoresCombatOrigDupSeed$impute, drop = TRUE)
# 
# scoresGender$method <- interaction(scoresGender$cluster,
#                                    scoresGender$impute, drop = TRUE)
# 
# scoresGenderOrig$method <- interaction(scoresGenderOrig$cluster,
#                                        scoresGenderOrig$impute, drop = TRUE)
# 
# # separate into male and female for complete data
# scoresGenderMale <- scoresGender[scoresGender$gender == 1,]
# scoresGenderFemale <- scoresGender[scoresGender$gender == 2,]
# 
# # separate into male and female
# scoresGenderOrigMale <- scoresGenderOrig[scoresGenderOrig$gender == 1,]
# scoresGenderOrigFemale <- scoresGenderOrig[scoresGenderOrig$gender == 2,]
# 
# scoresNormal3000$method <- interaction(scoresNormal3000$cluster,
#                                       scoresNormal3000$impute, drop = TRUE)
# 
# scoresUnion3000$method <- interaction(scoresUnion3000$cluster,
#                                        scoresUnion3000$impute, drop = TRUE)
# 

###################################################################################################
# Drop acc and nmi from original data, as it's not needed for a evaluation
scoresNormalOrig$acc <- NULL
scoresNormalOrig$nmi <- NULL
# scoresNormalOrigClust$acc <- NULL
# scoresNormalOrigClust$nmi <- NULL
# scoresCombatOrig$acc <- NULL
# scoresCombatOrig$nmi <- NULL
# scoresGenderOrig$acc <- NULL
# scoresGenderOrig$nmi <- NULL
scoresNormalOrigDup$acc <- NULL
scoresNormalOrigDup1000$nmi <- NULL
scoresNormalOrigDup1000$acc <- NULL
scoresNormalOrigDup3000$nmi <- NULL
scoresNormalOrigDup3000$acc <- NULL
scoresNormalOrigDup$nmi <- NULL
# scoresNormalOrigNA$acc <- NULL
# scoresNormalOrigNA$nmi <- NULL
# scoresUnion3000$acc <- NULL
# scoresUnion3000$nmi <- NULL

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
                          title) {
  
   if ((orig) && (pval)) {
      temp <- data %>%
        group_by(method) %>%
        summarise(meanPval = mean(pval, na.rm = T),
                  meanCi = mean(ci, na.rm = T))
      
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
# groupbyMethod(scoresNormalOrigDup, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union All Cancers Dup')
# groupbyMethod(scoresNormalOrigDup, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion All Cancers Dup')
# groupbyMethod(scoresNormalOrigDup, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union All Cancers Dup')
# groupbyMethod(scoresNormalOrigDup, orig = TRUE, bias = TRUE, title = 'Bias Union All Cancers Dup')
# groupbyMethod(scoresNormalOrigDup, orig = TRUE, coef = TRUE, title = 'Coef Union All Cancers Dup')

groupbyMethod(scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup 1000')
groupbyMethod(scoresNormalOrigDup1000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union All Cancers Dup 1000')
# groupbyMethod(scoresNormalOrigDup1000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion All Cancers Dup 1000')
# groupbyMethod(scoresNormalOrigDup1000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union All Cancers Dup 1000')
# groupbyMethod(scoresNormalOrigDup1000, orig = TRUE, bias = TRUE, title = 'Bias Union All Cancers Dup 1000')
# groupbyMethod(scoresNormalOrigDup1000, orig = TRUE, coef = TRUE, title = 'Coef Union All Cancers Dup 1000')

groupbyMethod(scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers Dup 3000')
# groupbyMethod(scoresNormalOrigDup3000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union All Cancers Dup 3000')
# groupbyMethod(scoresNormalOrigDup3000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion All Cancers Dup 3000')
# groupbyMethod(scoresNormalOrigDup3000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union All Cancers Dup 3000')
# groupbyMethod(scoresNormalOrigDup3000, orig = TRUE, bias = TRUE, title = 'Bias Union All Cancers Dup 3000')
# groupbyMethod(scoresNormalOrigDup3000, orig = TRUE, coef = TRUE, title = 'Coef Union All Cancers Dup 3000')

# groupbyMethod(scoresNormalOrigDupSeed, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers DupSeed')
# groupbyMethod(scoresNormalOrigDupSeed, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union All Cancers DupSeed')
# groupbyMethod(scoresNormalOrigDupSeed, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion All Cancers DupSeed')
# groupbyMethod(scoresNormalOrigDupSeed, orig = TRUE, con_index_ci = TRUE, title = 'Con CiUnion All Cancers DupSeed')
# groupbyMethod(scoresNormalOrigDupSeed, orig = TRUE, bias = TRUE, title = 'Bias Union All Cancers DupSeed')
# groupbyMethod(scoresNormalOrigDupSeed, orig = TRUE, coef = TRUE, title = 'Coef Union All Cancers Dupseed')
# 
# groupbyMethod(scoresNormalOrigNA, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers NA')
# groupbyMethod(scoresNormalOrigNA, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union All Cancers NA')
# groupbyMethod(scoresNormalOrigNA, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion All Cancers NA')
# groupbyMethod(scoresNormalOrigNA, orig = TRUE, con_index_ci = TRUE, title = 'Con CiUnion All Cancers NA')
# groupbyMethod(scoresNormalOrigNA, orig = TRUE, bias = TRUE, title = 'Bias Union All Cancers NA')
# groupbyMethod(scoresNormalOrigNA, orig = TRUE, coef = TRUE, title = 'Coef Union All Cancers NA')
# 
# groupbyMethod(scoresNormalOrigClust, orig = TRUE, pval = TRUE, title = 'Original Data')
# 
# groupbyMethod(scoresGenderMale,normal = TRUE, title = 'Complete Data Male')
# groupbyMethod(scoresGenderOrigMale, orig = TRUE, normal = TRUE, title = 'Original Data Male')
# 
# groupbyMethod(scoresGenderFemale, normal = TRUE, title = 'Complete Data with Female')
# groupbyMethod(scoresGenderOrigFemale, orig = TRUE, pval = TRUE, title = 'Original Data with Female')
# 



# ############################################################################################
# # Group by cluster
# 
# groupbyCluster <- function(data, orig = FALSE, title) {
#   
#   if (orig) {
#     temp <- data %>%
#       group_by(cluster) %>%
#       summarise(meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     
#     temp_melt <- melt(temp, id.vars = c('cluster'))
#     ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       xlab('cluster') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
#     
#   } else {
#     
#     temp <- data %>%
#       group_by(cluster) %>%
#       summarise(meanAcc = mean(acc, na.rm = T),
#                 meanNmi = mean(nmi, na.rm = T),
#                 meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     
#     temp_melt <- melt(temp, id.vars = c('cluster'))
#     
#     ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       xlab('cluster') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
#   }
# }
# 
# 
# groupbyCluster(scoresNormal, title = 'Intersection All Cancers')
# groupbyCluster(scoresNormalOrig, orig = TRUE, title = 'Union All Cancers')
# groupbyCluster(scoresNormalOrigDupSeed, orig = TRUE, title = 'Union All Cancers with Seeds and duplications removed')
# groupbyCluster(scoresNormalOrigDup, orig = TRUE, title = 'Union No Duplicates All Cancers')
# groupbyCluster(scoresNormalOrigNA, orig = TRUE, title = 'Original Data')
# groupbyCluster(scoresNormalOrigClust, orig = TRUE, title = 'Original Data')
# groupbyCluster(scoresLUSCOrigDup, orig = TRUE, title = 'LUSC Original Data')
# 
# groupbyCluster(scoresCombat, title = 'Complete Data with Combat')
# groupbyCluster(scoresCombatOrig, orig = TRUE, title = 'Original Data with Combat')
# groupbyCluster(scoresCombatOrigDup, orig = TRUE, title = 'Original Data with Combat, Revoed Dups')
# groupbyCluster(scoresCombatOrigDupSeed, orig = TRUE, title = 'Original Data with Combat, Removed Dups and random')
# 
# groupbyCluster(scoresGenderMale, title = 'Complete Data Male')
# groupbyCluster(scoresGenderOrigMale, orig = TRUE, title = 'Original Data Male')
# 
# groupbyCluster(scoresGenderFemale, title = 'Complete Data with Female')
# groupbyCluster(scoresGenderOrigFemale, orig = TRUE, title = 'Original Data with Female')
# 
# ############################################################################################
# # Group by impute
# 
# groupbyimpute <- function(data, orig = FALSE, title) {
#   
#   if (orig) {
#     temp <- data %>%
#       group_by(impute) %>%
#       summarise(meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     
#     temp_melt <- melt(temp, id.vars = c('impute'))
#     ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       xlab('impute') + ylab('Mean Score') +  ggtitle(title) + theme_538_bar
#     
#   } else {
#     
#     temp <- data %>%
#       group_by(impute) %>%
#       summarise(meanAcc = mean(acc, na.rm = T),
#                 meanNmi = mean(nmi, na.rm = T),
#                 meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     
#     temp_melt <- melt(temp, id.vars = c('impute'))
#     
#     ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       xlab('impute') + ylab('Mean Score') + ggtitle(title) + theme_538_bar
#   }
# }
# 
# 
# groupbyimpute(scoresNormal, title = 'Intersection All Cancers')
# groupbyimpute(scoresNormalOrig, orig = TRUE, title = 'Union All Cancers')
# groupbyimpute(scoresNormalOrigDupSeed, orig = TRUE, title = 'Union All Cancers with Seeds and duplications removed')
# groupbyimpute(scoresNormalOrigDup, orig = TRUE, title = 'Union No Duplicates All Cancers')
# groupbyimpute(scoresNormalOrigNA, orig = TRUE, title = 'Original Data')
# groupbyimpute(scoresNormalOrigClust, orig = TRUE, title = 'Original Data')
# groupbyimpute(scoresLUSCOrigDup, orig = TRUE, title = 'LUSC Original Data')
# 
# 
# groupbyimpute(scoresCombat, title = 'Complete Data with Combat')
# groupbyimpute(scoresCombatOrig, orig = TRUE, title = 'Original Data with Combat')
# groupbyimpute(scoresCombatOrigDup, orig = TRUE, title = 'Original Data with Combat, Revoed Dups')
# groupbyimpute(scoresCombatOrigDupSeed, orig = TRUE, title = 'Original Data with Combat, Removed Dups and random')
# 
# 
# 
# groupbyimpute(scoresGenderMale, title = 'Complete Data Male')
# groupbyimpute(scoresGenderOrigMale, orig = TRUE, title = 'Original Data Male')
# 
# groupbyimpute(scoresGenderFemale, title = 'Complete Data with Female')
# groupbyimpute(scoresGenderOrigFemale, orig = TRUE, title = 'Original Data with Female')

###############################
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
  }
}

groupbyCancer(cancer = 1, scoresNormal, normal = TRUE, title = 'Intersection BRCA')

groupbyCancer(cancer = 1, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA')

groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA Dup')
# groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union BRCA Dup')
# groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion BRCA Dup')
# groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union BRCA Dup')
# groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, bias = TRUE, title = 'bias Union BRCA Dup')
# groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, coef = TRUE, title = 'coef Union BRCA Dup')

groupbyCancer(cancer = 1, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA Dup 1000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup1000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union BRCA Dup 1000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup1000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion BRCA Dup 1000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup1000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union BRCA Dup 1000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup1000, orig = TRUE, bias = TRUE, title = 'bias Union BRCA Dup 1000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup1000, orig = TRUE, coef = TRUE, title = 'coef Union BRCA Dup 1000')

groupbyCancer(cancer = 1, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA Dup 3000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup3000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union BRCA Dup 3000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup3000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion BRCA Dup 3000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup3000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union BRCA Dup 3000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup3000, orig = TRUE, bias = TRUE, title = 'bias Union BRCA Dup 3000')
# groupbyCancer(cancer = 1, scoresNormalOrigDup3000, orig = TRUE, coef = TRUE, title = 'coef Union BRCA Dup 3000')

# groupbyCancer(cancer = 1, scoresNormalOrigDupSeed, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA DupSeed')
# groupbyCancer(cancer = 1, scoresNormalOrigDupSeed, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union BRCA DupSeed')
# groupbyCancer(cancer = 1, scoresNormalOrigDupSeed, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion BRCA DupSeed')
# groupbyCancer(cancer = 1, scoresNormalOrigDupSeed, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union BRCA DupSeed')
# groupbyCancer(cancer = 1, scoresNormalOrigDupSeed, orig = TRUE, bias = TRUE, title = 'bias Union BRCA DupSeed')
# groupbyCancer(cancer = 1, scoresNormalOrigDupSeed, orig = TRUE, coef = TRUE, title = 'coef Union BRCA DupSeed')
# 
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, pval = TRUE, title = 'Pval Union BRCA NA')
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union BRCA NA')
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion BRCA Na')
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union BRCA NA')
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, bias = TRUE, title = 'bias Union BRCA NA')
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, coef = TRUE, title = 'coef Union BRCA NA')

## KIRC

groupbyCancer(cancer = 2, scoresNormal, normal = TRUE, title = 'Intersection KIRC')

groupbyCancer(cancer = 2, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union kirc')

groupbyCancer(cancer = 2, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union kirc Dup')

groupbyCancer(cancer = 2, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union kirc Dup 1000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup1000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union kirc Dup 1000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup1000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion kirc Dup 1000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup1000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union kirc Dup 1000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup1000, orig = TRUE, bias = TRUE, title = 'bias Union kirc Dup 1000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup1000, orig = TRUE, coef = TRUE, title = 'coef Union kirc Dup 1000')

# groupbyCancer(cancer = 2, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union kirc Dup 3000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup3000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union kirc Dup 3000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup3000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion kirc Dup 3000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup3000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union kirc Dup 3000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup3000, orig = TRUE, bias = TRUE, title = 'bias Union kirc Dup 3000')
# groupbyCancer(cancer = 2, scoresNormalOrigDup3000, orig = TRUE, coef = TRUE, title = 'coef Union kirc Dup 3000')

# groupbyCancer(cancer = 2, scoresNormalOrigDupSeed, orig = TRUE, pval = TRUE, title = 'Pval Union kirc DupSeed')
# groupbyCancer(cancer = 2, scoresNormalOrigDupSeed, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union kirc DupSeed')
# groupbyCancer(cancer = 2, scoresNormalOrigDupSeed, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion kirc DupSeed')
# groupbyCancer(cancer = 2, scoresNormalOrigDupSeed, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union kirc DupSeed')
# groupbyCancer(cancer = 2, scoresNormalOrigDupSeed, orig = TRUE, bias = TRUE, title = 'bias Union kirc DupSeed')
# groupbyCancer(cancer = 2, scoresNormalOrigDupSeed, orig = TRUE, coef = TRUE, title = 'coef Union kirc DupSeed')
# 
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, pval = TRUE, title = 'Pval Union kirc NA')
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union kirc NA')
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion kirc Na')
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union kirc NA')
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, bias = TRUE, title = 'bias Union kirc NA')
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, coef = TRUE, title = 'coef Union kirc NA')


groupbyMethod(scoresCombat, normal = TRUE, title = 'Combat intersection')
groupbyMethod(scoresCombatDup, normal = TRUE, title = 'Combat intersection')


# groupbyMethod(scoresCombatOrig, orig = TRUE, pval = TRUE, title = 'Combat union')

groupbyMethod(scoresCombatOrigDup, orig = TRUE, pval = TRUE, title = 'Combat union duplicate removed')
# groupbyMethod(scoresCombatOrigDup, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union combat DupSeed')
# groupbyMethod(scoresCombatOrigDup, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion combat DupSeed')
# groupbyMethod(scoresCombatOrigDup, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union combat DupSeed')
# groupbyMethod(scoresCombatOrigDup, orig = TRUE, bias = TRUE, title = 'bias Union combat DupSeed')
# groupbyMethod(scoresCombatOrigDup, orig = TRUE, coef = TRUE, title = 'coef Union combat DupSeed')

# groupbyMethod(scoresCombatOrigDupSeed, orig = TRUE, pval = TRUE, title = 'Pval Union combat DupSeed')
# groupbyMethod(scoresCombatOrigDupSeed, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union combat DupSeed')
# groupbyMethod(scoresCombatOrigDupSeed, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion combat DupSeed')
# groupbyMethod(scoresCombatOrigDupSeed, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union combat DupSeed')
# groupbyMethod(scoresCombatOrigDupSeed, orig = TRUE, bias = TRUE, title = 'bias Union combat DupSeed')
# groupbyMethod(scoresCombatOrigDupSeed, orig = TRUE, coef = TRUE, title = 'coef Union combat DupSeed')

## lihc
groupbyCancer(cancer = 3, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union lihc')

groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union lihc Dup')
# groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union lihc Dup')
# groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion lihc Dup')
# groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union lihc Dup')
# groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, bias = TRUE, title = 'bias Union lihc Dup')
# groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, coef = TRUE, title = 'coef Union lihc Dup')

groupbyCancer(cancer = 3, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union lihc Dup 1000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup1000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union lihc Dup 1000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup1000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion lihc Dup 1000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup1000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union lihc Dup 1000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup1000, orig = TRUE, bias = TRUE, title = 'bias Union lihc Dup 1000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup1000, orig = TRUE, coef = TRUE, title = 'coef Union lihc Dup 1000')

# groupbyCancer(cancer = 3, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union lihc Dup 3000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup3000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union lihc Dup 3000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup3000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion lihc Dup 3000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup3000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union lihc Dup 3000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup3000, orig = TRUE, bias = TRUE, title = 'bias Union lihc Dup 3000')
# groupbyCancer(cancer = 3, scoresNormalOrigDup3000, orig = TRUE, coef = TRUE, title = 'coef Union lihc Dup 3000')

# groupbyCancer(cancer = 3, scoresNormalOrigDupSeed, orig = TRUE, pval = TRUE, title = 'Pval Union lihc DupSeed')
# groupbyCancer(cancer = 3, scoresNormalOrigDupSeed, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union lihc DupSeed')
# groupbyCancer(cancer = 3, scoresNormalOrigDupSeed, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion lihc DupSeed')
# groupbyCancer(cancer = 3, scoresNormalOrigDupSeed, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union lihc DupSeed')
# groupbyCancer(cancer = 3, scoresNormalOrigDupSeed, orig = TRUE, bias = TRUE, title = 'bias Union lihc DupSeed')
# groupbyCancer(cancer = 3, scoresNormalOrigDupSeed, orig = TRUE, coef = TRUE, title = 'coef Union lihc DupSeed')
# 
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, pval = TRUE, title = 'Pval Union lihc NA')
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union lihc NA')
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion lihc Na')
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union lihc NA')
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, bias = TRUE, title = 'bias Union lihc NA')
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, coef = TRUE, title = 'coef Union lihc NA')
# 

## luad
groupbyCancer(cancer = 4, scoresNormalOrig, orig = TRUE, pval = TRUE, title = 'Pval Union luad')

groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, pval = TRUE, title = 'Pval Union luad Dup')
groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union luad Dup')
groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion luad Dup')
groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union luad Dup')
groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, bias = TRUE, title = 'bias Union luad Dup')
groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, coef = TRUE, title = 'coef Union luad Dup')

groupbyCancer(cancer = 4, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union luad Dup 1000')
groupbyCancer(cancer = 4, scoresNormalOrigDup1000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union luad Dup 1000')
groupbyCancer(cancer = 4, scoresNormalOrigDup1000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion luad Dup 1000')
groupbyCancer(cancer = 4, scoresNormalOrigDup1000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union luad Dup 1000')
groupbyCancer(cancer = 4, scoresNormalOrigDup1000, orig = TRUE, bias = TRUE, title = 'bias Union luad Dup 1000')
groupbyCancer(cancer = 4, scoresNormalOrigDup1000, orig = TRUE, coef = TRUE, title = 'coef Union luad Dup 1000')

groupbyCancer(cancer = 4, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union luad Dup 3000')
groupbyCancer(cancer = 4, scoresNormalOrigDup3000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union luad Dup 3000')
groupbyCancer(cancer = 4, scoresNormalOrigDup3000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion luad Dup 3000')
groupbyCancer(cancer = 4, scoresNormalOrigDup3000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union luad Dup 3000')
groupbyCancer(cancer = 4, scoresNormalOrigDup3000, orig = TRUE, bias = TRUE, title = 'bias Union luad Dup 3000')
groupbyCancer(cancer = 4, scoresNormalOrigDup3000, orig = TRUE, coef = TRUE, title = 'coef Union luad Dup 3000')

groupbyCancer(cancer = 4, scoresNormalOrigDupSeed, orig = TRUE, pval = TRUE, title = 'Pval Union luad DupSeed')
groupbyCancer(cancer = 4, scoresNormalOrigDupSeed, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union luad DupSeed')
groupbyCancer(cancer = 4, scoresNormalOrigDupSeed, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion luad DupSeed')
groupbyCancer(cancer = 4, scoresNormalOrigDupSeed, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union luad DupSeed')
groupbyCancer(cancer = 4, scoresNormalOrigDupSeed, orig = TRUE, bias = TRUE, title = 'bias Union luad DupSeed')
groupbyCancer(cancer = 4, scoresNormalOrigDupSeed, orig = TRUE, coef = TRUE, title = 'coef Union luad DupSeed')

groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, pval = TRUE, title = 'Pval Union luad NA')
groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union luad NA')
groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion luad Na')
groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union luad NA')
groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, bias = TRUE, title = 'bias Union luad NA')
groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, coef = TRUE, title = 'coef Union luad NA')

# LUSC

groupbyCancer(cancer = 5, scoresNormalOrigDup1000, orig = TRUE, pval = TRUE, title = 'Pval Union lusc Dup 1000')
groupbyCancer(cancer = 5, scoresNormalOrigDup1000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union lusc Dup 1000')
groupbyCancer(cancer = 5, scoresNormalOrigDup1000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion lusc Dup 1000')
groupbyCancer(cancer = 5, scoresNormalOrigDup1000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union lusc Dup 1000')
groupbyCancer(cancer = 5, scoresNormalOrigDup1000, orig = TRUE, bias = TRUE, title = 'bias Union lusc Dup 1000')
groupbyCancer(cancer = 5, scoresNormalOrigDup1000, orig = TRUE, coef = TRUE, title = 'coef Union lusc Dup 1000')

groupbyCancer(cancer = 5, scoresNormalOrigDup3000, orig = TRUE, pval = TRUE, title = 'Pval Union lusc Dup 3000')
groupbyCancer(cancer = 5, scoresNormalOrigDup3000, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union lusc Dup 3000')
groupbyCancer(cancer = 5, scoresNormalOrigDup3000, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion lusc Dup 3000')
groupbyCancer(cancer = 5, scoresNormalOrigDup3000, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union lusc Dup 3000')
groupbyCancer(cancer = 5, scoresNormalOrigDup3000, orig = TRUE, bias = TRUE, title = 'bias Union lusc Dup 3000')
groupbyCancer(cancer = 5, scoresNormalOrigDup3000, orig = TRUE, coef = TRUE, title = 'coef Union lusc Dup 3000')


groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, pval = TRUE, title = 'Intersection LUSC')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, pvalcox = TRUE, title = 'Pvalcox Union luad Dup')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, con_index_p = TRUE, title = 'Con PUnion luad Dup')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, con_index_ci = TRUE, title = 'Con Ci Union luad Dup')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, bias = TRUE, title = 'bias Union luad Dup')
groupbyCancer(cancer = 5, scoresLUSCNormalDup, normal = TRUE, coef = TRUE, title = 'coef Union luad Dup')

groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, pval = TRUE, title = 'Union LUSC duplicates removed')
groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, pvalcox = TRUE, title = 'Pvalcox Union luad Dup')
groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, con_index_p = TRUE, title = 'Con PUnion luad Dup')
groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, con_index_ci = TRUE, title = 'Con Ci Union luad Dup')
groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, bias = TRUE, title = 'bias Union luad Dup')
groupbyCancer(cancer = 5, scoresLUSCOrigDup, orig = TRUE, coef = TRUE, title = 'coef Union luad Dup')



# ##################################################################################################
# # Group by cancer and method and mean of evaluation. 
# groupbyCancer <- function(cancer, data, orig = FALSE, title) {
#   
#   temp_data <- data[data$cancer == cancer,]
#   
#   if (orig) {
#     
#     temp <- temp_data %>%
#       group_by(cluster) %>%
#       summarise(meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     temp_melt <- melt(temp, id.vars = c('cluster'))
#     ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('cluster') + ggtitle(title)
#     
#   } else {
#     
#     temp <- temp_data %>%
#       group_by(cluster) %>%
#       summarise(meanAcc = mean(acc, na.rm = T),
#                 meanNmi = mean(nmi, na.rm = T),
#                 meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     temp_melt <- melt(temp, id.vars = c('cluster'))
#     
#     ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('cluster') + ggtitle(title)
#   }
# }
# 
# groupbyCancer(cancer = 1, scoresNormal, title = 'BRCA Complete Data')
# groupbyCancer(cancer = 2, scoresNormal, title = 'KIRC Complete Data')
# groupbyCancer(cancer = 3, scoresNormal, title = 'LIHC Complete Data')
# groupbyCancer(cancer = 4, scoresNormal, title = 'LUAD Complete Data')
# 
# groupbyCancer(cancer = 1, scoresNormalOrig, orig = TRUE, title = 'BRCA Original Data')
# groupbyCancer(cancer = 2, scoresNormalOrig, orig = TRUE, title = 'KIRC Original Data')
# groupbyCancer(cancer = 3, scoresNormalOrig, orig = TRUE, title = 'LIHC Original Data')
# groupbyCancer(cancer = 4, scoresNormalOrig, orig = TRUE, title = 'LUAD Original Data')
# 
# groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, title = 'BRCA Original Data')
# groupbyCancer(cancer = 2, scoresNormalOrigDup, orig = TRUE, title = 'KIRC Original Data')
# groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, title = 'LIHC Original Data')
# groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, title = 'LUAD Original Data')
# 
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, title = 'BRCA Original Data')
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, title = 'KIRC Original Data')
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, title = 'LIHC Original Data')
# groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, title = 'LUAD Original Data')
# 
# 
# 
# ##################################################################################################
# # Group by cancer and method and mean of evaluation. 
# groupbyCancer <- function(cancer, data, orig = FALSE, title) {
#   
#   temp_data <- data[data$cancer == cancer,]
#   
#   if (orig){
#     
#     temp <- temp_data %>%
#       group_by(impute) %>%
#       summarise(meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     temp_melt <- melt(temp, id.vars = c('impute'))
#     ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('impute') + ggtitle(title)
#     
#   } else {
#     
#     temp <- temp_data %>%
#       group_by(impute) %>%
#       summarise(meanAcc = mean(acc, na.rm = T),
#                 meanNmi = mean(nmi, na.rm = T),
#                 meanPval = mean(pval, na.rm = T),
#                 meanCi = mean(ci, na.rm = T))
#     temp_melt <- melt(temp, id.vars = c('impute'))
#     
#     ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('impute') + ggtitle(title)
#   }
# }
# 
# groupbyCancer(cancer = 1, scoresNormal, title = 'BRCA Complete Data')
# groupbyCancer(cancer = 2, scoresNormal, title = 'KIRC Complete Data')
# groupbyCancer(cancer = 3, scoresNormal, title = 'LIHC Complete Data')
# groupbyCancer(cancer = 4, scoresNormal, title = 'LUAD Complete Data')
# 
# groupbyCancer(cancer = 1, scoresNormalOrig, orig = TRUE, title = 'BRCA Original Data')
# groupbyCancer(cancer = 2, scoresNormalOrig, orig = TRUE, title = 'KIRC Original Data')
# groupbyCancer(cancer = 3, scoresNormalOrig, orig = TRUE, title = 'LIHC Original Data')
# groupbyCancer(cancer = 4, scoresNormalOrig, orig = TRUE, title = 'LUAD Original Data')
# 
# groupbyCancer(cancer = 1, scoresNormalOrigDup, orig = TRUE, title = 'BRCA Original Data')
# groupbyCancer(cancer = 2, scoresNormalOrigDup, orig = TRUE, title = 'KIRC Original Data')
# groupbyCancer(cancer = 3, scoresNormalOrigDup, orig = TRUE, title = 'LIHC Original Data')
# groupbyCancer(cancer = 4, scoresNormalOrigDup, orig = TRUE, title = 'LUAD Original Data')
# 
# groupbyCancer(cancer = 1, scoresNormalOrigNA, orig = TRUE, title = 'BRCA Original Data')
# groupbyCancer(cancer = 2, scoresNormalOrigNA, orig = TRUE, title = 'KIRC Original Data')
# groupbyCancer(cancer = 3, scoresNormalOrigNA, orig = TRUE, title = 'LIHC Original Data')
# groupbyCancer(cancer = 4, scoresNormalOrigNA, orig = TRUE, title = 'LUAD Original Data')
# 
# 

###############################################################################################################
# Distribution of acc, nmi, pval, and ci
distComplete <- function(column) {
  par(mfrow = c(2,2))
  hist(scoresNormal[, column], main = paste0('complete normal', ' ', column), 
       xlab = 'Normal', col = 'lightblue')
  hist(scoresCombat[, column], main = paste0('complete combat', ' ', column), 
       xlab = 'Combat', col = 'lightblue' )
  hist(scoresGenderMale[, column], main = paste0('complete male', ' ', column),
       xlab = 'Male', col = 'lightblue')
  hist(scoresGenderFemale[, column], main = paste0('complete female', ' ', column), 
       xlab = 'Female', col = 'lightblue')
}
distComplete('acc')
distComplete('nmi')
distComplete('pval')
distComplete('ci')


distOrig <- function(column) {
  par(mfrow =  c(2,2))
  hist(scoresNormalOrig[, column], main = paste0('original normal', ' ', column), xlab = 'Normal', col = 'lightblue')
  hist(scoresCombatOrig[, column], main = paste0('original combat', ' ', column), xlab = 'Combat', col = 'lightblue')
  hist(scoresGenderOrigMale[, column], main = paste0('original male', ' ', column), xlab = 'Male', col = 'lightblue')
  hist(scoresGenderOrigFemale[, column], main = paste0('original female', ' ', column), xlab = 'Female', col = 'lightblue')
}
distOrig('pval')
distOrig('ci')

#########################################################################################################
# Actual pval 

groupbyPval <- function(data, data2, data_indicator = 'Normal', cancer = NULL, title) {
  
  if (grepl('Combat', data_indicator)) {
    temp <- data %>%
      group_by(method) %>%
      summarise(actual_pval_combat = mean(actual_pval, na.rm = T))
    data2 <- data2[which(data2$cancer == 2),]
   
    temp.2 <- data2 %>%
      group_by(method) %>%
      summarise(actual_pval_normal = mean(actual_pval, na.rm = T))
    temp.full <- merge(temp, temp.2, by = 'method')
    
    temp_melt <- melt(temp.full, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, value), value, fill = variable)) + geom_bar(stat = 'identity') +
     xlab('Method') + ggtitle(title) + theme_538_bar
  
    } else {
  
    data <- data[which(data$cancer == cancer),]
    temp <- data %>%
      group_by(method) %>%
      summarise(actual_pval = mean(actual_pval, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, value), value, fill = variable)) + geom_bar(stat = 'identity') +
      xlab('Method') + ggtitle(title) + theme_538_bar
  }

}

groupbyPval(scoresNormal, title = 'Complete Data')
groupbyPval(scoresNormalOrig, title = 'Union Data')



groupbyPval(scoresNormal, title = 'Complete Data BRCA', cancer = 1)
groupbyPval(scoresNormal, title = 'Complete Data KIRC', cancer = 2)
groupbyPval(scoresNormal, title = 'Complete Data LIHC', cancer = 3)
groupbyPval(scoresNormal, title = 'Complete Data LUAD', cancer = 4)

groupbyPval(scoresNormalOrig, title = 'Union Data BRCA', cancer = 1)
groupbyPval(scoresNormalOrig, title = 'Union Data KIRC', cancer = 2)
groupbyPval(scoresNormalOrig, title = 'Union Data LIHC', cancer = 3)
groupbyPval(scoresNormalOrig, title = 'Union Data LUAD', cancer = 4)

groupbyPval(scoresNormalOrigDupSeed, title = 'Union Data BRCA Dup/Seed', cancer = 1)
groupbyPval(scoresNormalOrigDupSeed, title = 'Union Data KIRC Dup/Seed', cancer = 2)
groupbyPval(scoresNormalOrigDupSeed, title = 'Union Data LIHC Dup/Seed' cancer = 3)
groupbyPval(scoresNormalOrigDupSeed, title = 'Union Data LUAD Dup/Seed', cancer = 4)
groupbyPval(scoresLUSCOrigDup, title = 'Union Data LUSC Dup/Seed', cancer = 5)


groupbyPval(scoresNormalOrigClust, title = 'Union Data BRCA Clust', cancer = 1)
groupbyPval(scoresNormalOrigClust, title = 'Union Data KIRC Clust', cancer = 2)
groupbyPval(scoresNormalOrigClust, title = 'Union Data LIHC Clust', cancer = 3)
groupbyPval(scoresNormalOrigClust, title = 'Union Data LUAD Clust', cancer = 4)


groupbyPval(scoresNormalOrig, title = 'Original Data')
groupbyPval(scoresNormalOrigClust, title = 'Original Data')


groupbyPval(scoresCombat, title = 'KIRC Complete with Combat')
groupbyPval(scoresCombatOrig, title = 'KIRC Original with Combat')
groupbyPval(scoresCombatOrigDupSeed, title = 'KIRC Complete with Combat/Seed')
groupbyPval(scoresCombatOrig, title = 'KIRC Original with Combat')
groupbyPval(scoresCombat, scoresNormal, data_indicator = 'Combat', title = 'KIRC Complete with Combat and Normal')
groupbyPval(scoresCombatOrig, scoresNormalOrig, data_indicator = 'Combat', title = 'KIRC Original with Combat and Normal')

groupbyPval(scoresGenderMale, title = 'Complete Data Male')
groupbyPval(scoresGenderOrigMale, title = 'Original Data Male')

groupbyPval(scoresGenderFemale, title = 'Complete Data with Female')
groupbyPval(scoresGenderOrigFemale,title = 'Original Data with Female')

################################################################################################################
# Look at intersection of union
temp <- scoresNormalOrigInt %>%
  group_by(method) %>%
  summarise(meanAcc = mean(intAcc, na.rm = T),
            meanNmi = mean(intNmi, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('Intersection of Union') + theme_538_bar

#Look at intersection
temp <- scoresNormalOrigInt %>%
  group_by(cluster) %>%
  summarise(meanAcc = mean(intAcc, na.rm = T),
            meanNmi = mean(intNmi, na.rm = T))

temp_melt <- melt(temp, id.vars = c('cluster'))

ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('Intersection of Union') + theme_538_bar

#Look at intersection
temp <- scoresNormalOrigInt %>%
  group_by(impute) %>%
  summarise(meanAcc = mean(intAcc, na.rm = T),
            meanNmi = mean(intNmi, na.rm = T))

temp_melt <- melt(temp, id.vars = c('impute'))

ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('Intersection of Union') + theme_538_bar


## By Cancer
groupByCancerInt <- function(cancer, title) {
  temp.data <- scoresNormalOrigInt[scoresNormalOrigInt$cancer == cancer,]
  
  temp <- temp.data %>%
    group_by(method) %>%
    summarise(meanAcc = mean(intAcc, na.rm = T),
              meanNmi = mean(intNmi, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('method'))
  
  ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title) + theme_538_bar
  
}

groupByCancerInt(1, title = 'BRCA Intersection of Union')
groupByCancerInt(2, title = 'KIRC Intersection of Union')
groupByCancerInt(3, title = 'LIHC Intersection of Union')
groupByCancerInt(4, title = 'LUAD Intersection of Union')


## By Cancer
groupByCancerInt <- function(cancer, title) {
  temp.data <- scoresNormalOrigInt[scoresNormalOrigInt$cancer == cancer,]
  
  temp <- temp.data %>%
    group_by(cluster) %>%
    summarise(meanAcc = mean(intAcc, na.rm = T),
              meanNmi = mean(intNmi, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('cluster'))
  
  ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title) + theme_538_bar
  
}

groupByCancerInt(1, title = 'BRCA Intersection of Union')
groupByCancerInt(2, title = 'KIRC Intersection of Union')
groupByCancerInt(3, title = 'LIHC Intersection of Union')
groupByCancerInt(4, title = 'LUAD Intersection of Union')


## By Cancer
groupByCancerInt <- function(cancer, title) {
  temp.data <- scoresNormalOrigInt[scoresNormalOrigInt$cancer == cancer,]
  
  temp <- temp.data %>%
    group_by(impute) %>%
    summarise(meanAcc = mean(intAcc, na.rm = T),
              meanNmi = mean(intNmi, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('impute'))
  
  ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title) + theme_538_bar
  
}

groupByCancerInt(1, title = 'BRCA Intersection of Union')
groupByCancerInt(2, title = 'KIRC Intersection of Union')
groupByCancerInt(3, title = 'LIHC Intersection of Union')
groupByCancerInt(4, title = 'LUAD Intersection of Union')


################################################################################################################
# Look at intersection of union
temp <- scoresNormalOrigIntDup %>%
  group_by(method) %>%
  summarise(meanAcc = mean(intAcc, na.rm = T),
            meanNmi = mean(intNmi, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('Intersection of Union') + theme_538_bar

#Look at intersection
temp <- scoresNormalOrigIntDup %>%
  group_by(cluster) %>%
  summarise(meanAcc = mean(intAcc, na.rm = T),
            meanNmi = mean(intNmi, na.rm = T))

temp_melt <- melt(temp, id.vars = c('cluster'))

ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('Intersection of Union') + theme_538_bar

#Look at intersection
temp <- scoresNormalOrigIntDup %>%
  group_by(impute) %>%
  summarise(meanAcc = mean(intAcc, na.rm = T),
            meanNmi = mean(intNmi, na.rm = T))

temp_melt <- melt(temp, id.vars = c('impute'))

ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('Intersection of Union') + theme_538_bar


## By Cancer
groupByCancerInt <- function(cancer, title) {
  temp.data <- scoresNormalOrigIntDup[scoresNormalOrigIntDup$cancer == cancer,]
  
  temp <- temp.data %>%
    group_by(method) %>%
    summarise(meanAcc = mean(intAcc, na.rm = T),
              meanNmi = mean(intNmi, na.rm = T))
  
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
    group_by(cluster) %>%
    summarise(meanAcc = mean(intAcc, na.rm = T),
              meanNmi = mean(intNmi, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('cluster'))
  
  ggplot(data = temp_melt, aes(reorder(cluster, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
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
    group_by(impute) %>%
    summarise(meanAcc = mean(intAcc, na.rm = T),
              meanNmi = mean(intNmi, na.rm = T))
  
  temp_melt <- melt(temp, id.vars = c('impute'))
  
  ggplot(data = temp_melt, aes(reorder(impute, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title) + theme_538_bar
  
}

groupByCancerInt(1, title = 'BRCA Intersection of Union')
groupByCancerInt(2, title = 'KIRC Intersection of Union')
groupByCancerInt(3, title = 'LIHC Intersection of Union')
groupByCancerInt(4, title = 'LUAD Intersection of Union')
groupByCancerInt(5, title = 'LUSC Intersection of Union')


######################################################################################################################
# different measure of pvalue and ci

# ########################################################################################
# complete data 
scoresNormalPval <- read.csv(paste0(results_folder, '/scoresTwoThousandPval.csv'))

scoresNormalPval$method <- interaction(scoresNormalPval$cluster,
                                       scoresNormalPval$impute, drop = TRUE)

# across cancer pvalcox
temp <- scoresNormalPval %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPvalcox = mean(pvalcox, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar


# across cancer con_index_p
temp <- scoresNormalPval %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanCon_index_p = mean(con_index_p, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer con_index_ci
temp <- scoresNormalPval %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCon_index_ci = mean(con_index_ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index_ci') + theme_538_bar


# across cancer bias_corrected_c_i
temp <- scoresNormalPval %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer coef
temp <- scoresNormalPval %>%
  group_by(method) %>%
  summarise(meancoef = mean(coefcox, na.rm = T),
            meanse.std = mean(se.std, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('coef') + theme_538_bar


########################################################################## BRCA 
brca <- scoresNormalPval[scoresNormalPval$cancer == 1,]
# across cancer pvalcox
temp <- brca %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPvalcox = mean(pvalcox, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar


# across cancer con_index_p
temp <- brca %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanCon_index_p = mean(con_index_p, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer con_index_ci
temp <- brca %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCon_index_ci = mean(con_index_ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index_ci') + theme_538_bar


# across cancer bias_corrected_c_i
temp <- brca %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer coef
temp <- brca %>%
  group_by(method) %>%
  summarise(meancoef = mean(coefcox, na.rm = T),
            meanse.std = mean(se.std, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('coef') + theme_538_bar


# kirc##########################################################################

kirc <- scoresNormalPval[scoresNormalPval$cancer == 2,]
# across cancer pvalcox
temp <- kirc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar

# across cancer pvalcox
temp <- kirc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPvalcox = mean(pvalcox, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar


# across cancer con_index_p
temp <- kirc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanCon_index_p = mean(con_index_p, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer con_index_ci
temp <- kirc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCon_index_ci = mean(con_index_ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index_ci') + theme_538_bar


# across cancer bias_corrected_c_i
temp <- kirc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer coef
temp <- kirc %>%
  group_by(method) %>%
  summarise(meancoef = mean(coefcox, na.rm = T),
            meanse.std = mean(se.std, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('coef') + theme_538_bar

#############################################################################################
# lihc
lihc <- scoresNormalPval[scoresNormalPval$cancer == 3,]
# across cancer pvalcox
temp <- lihc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar

# across cancer pvalcox
temp <- lihc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPvalcox = mean(pvalcox, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar


# across cancer con_index_p
temp <- lihc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanCon_index_p = mean(con_index_p, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer con_index_ci
temp <- lihc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCon_index_ci = mean(con_index_ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index_ci') + theme_538_bar


# across cancer bias_corrected_c_i
temp <- lihc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer coef
temp <- lihc %>%
  group_by(method) %>%
  summarise(meancoef = mean(coefcox, na.rm = T),
            meanse.std = mean(se.std, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('coef') + theme_538_bar

# luad #############################################################################################

luad<- scoresNormalPval[scoresNormalPval$cancer == 4,]
# across cancer pvalcox
temp <- luad %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar

# across cancer pvalcox
temp <- luad %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPvalcox = mean(pvalcox, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar


# across cancer con_index_p
temp <- luad %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanCon_index_p = mean(con_index_p, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer con_index_ci
temp <- luad %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCon_index_ci = mean(con_index_ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index_ci') + theme_538_bar


# across cancer bias_corrected_c_i
temp <- luad %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer coef
temp <- luad %>%
  group_by(method) %>%
  summarise(meancoef = mean(coefcox, na.rm = T),
            meanse.std = mean(se.std, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('coef') + theme_538_bar


# lusc #############################################################################################

lusc<- scoresLUSCNormalDup
# across cancer pvalcox
temp <- lusc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar

# across cancer pvalcox
temp <- lusc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPvalcox = mean(pvalcox, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('pvalcox') + theme_538_bar


# across cancer con_index_p
temp <- lusc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanCon_index_p = mean(con_index_p, na.rm = T),
            meanCi = mean(ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer con_index_ci
temp <- lusc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCon_index_ci = mean(con_index_ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index_ci') + theme_538_bar


# across cancer bias_corrected_c_i
temp <- lusc %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer coef
temp <- lusc %>%
  group_by(method) %>%
  summarise(meancoef = mean(coefcox, na.rm = T),
            meanse.std = mean(se.std, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('coef') + theme_538_bar

# ##########################################################################################################
# # Original data
# 
# scoresNormalOrigPval <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandPval.csv'))
# 
# scoresNormalOrigPval$method <- interaction(scoresNormalOrigPval$cluster,
#                                        scoresNormalOrigPval$impute, drop = TRUE)
# 
# # across cancer pvalcox
# temp <- scoresNormalOrigPval %>%
#   group_by(method) %>%
#   summarise(meanPvalcox = mean(pvalcox, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# 
# # across cancer con_index_p
# temp <- scoresNormalOrigPval %>%
#   group_by(method) %>%
#   summarise(meanCon_index_p = mean(con_index_p, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer con_index_ci
# temp <- scoresNormalOrigPval %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCon_index_ci = mean(con_index_ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index_ci') + theme_538_bar
# 
# 
# # across cancer bias_corrected_c_i
# temp <- scoresNormalOrigPval %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer coef
# temp <- scoresNormalOrigPval %>%
#   group_by(method) %>%
#   summarise(meancoef = mean(coefcox, na.rm = T),
#             meanse.std = mean(se.std, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('coef') + theme_538_bar
# 
# 
# ########################################################################## BRCA 
# brca <- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 1,]
# # across cancer pvalcox
# temp <- brca %>%
#   group_by(method) %>%
#   summarise(meanPvalcox = mean(pvalcox, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# 
# # across cancer con_index_p
# temp <- brca %>%
#   group_by(method) %>%
#   summarise(meanCon_index_p = mean(con_index_p, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer con_index_ci
# temp <- brca %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCon_index_ci = mean(con_index_ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index_ci') + theme_538_bar
# 
# 
# # across cancer bias_corrected_c_i
# temp <- brca %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer coef
# temp <- brca %>%
#   group_by(method) %>%
#   summarise(meancoef = mean(coefcox, na.rm = T),
#             meanse.std = mean(se.std, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('coef') + theme_538_bar
# 
# 
# # kirc##########################################################################
# 
# kirc <- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 2,]
# # across cancer pvalcox
# temp <- kirc %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# # across cancer pvalcox
# temp <- kirc %>%
#   group_by(method) %>%
#   summarise(meanPvalcox = mean(pvalcox, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# # across cancer con_index_p
# temp <- kirc %>%
#   group_by(method) %>%
#   summarise(meanCon_index_p = mean(con_index_p, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer con_index_ci
# temp <- kirc %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCon_index_ci = mean(con_index_ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index_ci') + theme_538_bar
# 
# 
# # across cancer bias_corrected_c_i
# temp <- kirc %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer coef
# temp <- kirc %>%
#   group_by(method) %>%
#   summarise(meancoef = mean(coefcox, na.rm = T),
#             meanse.std = mean(se.std, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('coef') + theme_538_bar
# 
# #############################################################################################
# # lihc
# lihc <- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 3,]
# # across cancer pvalcox
# temp <- lihc %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# # across cancer pvalcox
# temp <- lihc %>%
#   group_by(method) %>%
#   summarise(meanPvalcox = mean(pvalcox, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# 
# # across cancer con_index_p
# temp <- lihc %>%
#   group_by(method) %>%
#   summarise(meanCon_index_p = mean(con_index_p, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer con_index_ci
# temp <- lihc %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCon_index_ci = mean(con_index_ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index_ci') + theme_538_bar
# 
# 
# # across cancer bias_corrected_c_i
# temp <- lihc %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer coef
# temp <- lihc %>%
#   group_by(method) %>%
#   summarise(meancoef = mean(coefcox, na.rm = T),
#             meanse.std = mean(se.std, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('coef') + theme_538_bar
# 
# # luad #############################################################################################
# 
# luad<- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 4,]
# # across cancer pvalcox
# temp <- luad %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# # across cancer pvalcox
# temp <- luad %>%
#   group_by(method) %>%
#   summarise(meanPvalcox = mean(pvalcox, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('pvalcox') + theme_538_bar
# 
# 
# # across cancer con_index_p
# temp <- luad %>%
#   group_by(method) %>%
#   summarise(meanCon_index_p = mean(con_index_p, na.rm = T),
#             meanCi = mean(ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer con_index_ci
# temp <- luad %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanCon_index_ci = mean(con_index_ci, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index_ci') + theme_538_bar
# 
# 
# # across cancer bias_corrected_c_i
# temp <- luad %>%
#   group_by(method) %>%
#   summarise(meanPval = mean(pval, na.rm = T),
#             meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('con_index') + theme_538_bar
# 
# # across cancer coef
# temp <- luad %>%
#   group_by(method) %>%
#   summarise(meancoef = mean(coefcox, na.rm = T),
#             meanse.std = mean(se.std, na.rm = T))
# 
# temp_melt <- melt(temp, id.vars = c('method'))
# 
# ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
#   xlab('Method') + ggtitle('coef') + theme_538_bar


