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

# Load data
scoresNormal <- read.csv(paste0(results_folder, '/scoresTwoThousand.csv'))
scoresNormalOrig <- read.csv(paste0(results_folder, '/scoresTwoThousandOrig.csv'))

scoresCombat <- read.csv(paste0(results_folder, '/scoresCombat.csv'))
scoresCombatOrig <- read.csv(paste0(results_folder, '/scoresCombatOrig.csv'))

scoresGender <- read.csv(paste0(results_folder, '/scoresGender.csv'))
scoresGenderOrig <- read.csv(paste0(results_folder, '/scoresGenderOrig.csv'))

###################################################################################################
# Get method types
scoresNormal$method <- interaction(scoresNormal$cluster,
                                   scoresNormal$impute, drop = TRUE)
scoresNormalOrig$method <- interaction(scoresNormalOrig$cluster,
                                   scoresNormalOrig$impute, drop = TRUE)

scoresCombat$method <- interaction(scoresCombat$cluster,
                                   scoresCombat$impute, drop = TRUE)
scoresCombatOrig$method <- interaction(scoresCombatOrig$cluster,
                                       scoresCombatOrig$impute, drop = TRUE)

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




###################################################################################################
# Drop acc and nmi from original data, as it's not needed for a evaluation
scoresNormalOrig$acc <- NULL
scoresNormalOrig$nmi <- NULL
scoresCombatOrig$acc <- NULL
scoresCombatOrig$nmi <- NULL
scoresGenderOrig$acc <- NULL
scoresGenderOrig$nmi <- NULL

####################################################################################################
# Group by method and get mean for each evaluation method (acc, nmi, pval, ci)

groupbyMethod <- function(data, orig = FALSE, title) {
  
   if (orig){
      temp <- data %>%
        group_by(method) %>%
        summarise(meanPval = mean(pval, na.rm = T),
                  meanCi = mean(ci, na.rm = T))
      
      temp_melt <- melt(temp, id.vars = c('method'))
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title)
      
  } else {
    
    temp <- data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title)
  }
}


groupbyMethod(scoresNormal, title = 'Complete Data')
groupbyMethod(scoresNormalOrig, orig = TRUE, title = 'Original Data')

groupbyMethod(scoresCombat, title = 'Complete Data with Combat')
groupbyMethod(scoresCombatOrig, orig = TRUE, title = 'Original Data with Combat')

groupbyMethod(scoresGenderMale, title = 'Complete Data Male')
groupbyMethod(scoresGenderOrigMale, orig = TRUE, title = 'Original Data Male')

groupbyMethod(scoresGenderFemale, title = 'Complete Data with Female')
groupbyMethod(scoresGenderOrigFemale, orig = TRUE, title = 'Original Data with Female')

##################################################################################################
# Group by cancer and method and mean of evaluation. 


groupbyCancer <- function(cancer, data, orig = FALSE, title) {
  
  temp_data <- data[data$cancer == cancer,]
  
  if (orig){
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    temp_melt <- melt(temp, id.vars = c('method'))
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title)
    
  } else {
    
    temp <- temp_data %>%
      group_by(method) %>%
      summarise(meanAcc = mean(acc, na.rm = T),
                meanNmi = mean(nmi, na.rm = T),
                meanPval = mean(pval, na.rm = T),
                meanCi = mean(ci, na.rm = T))
    temp_melt <- melt(temp, id.vars = c('method'))
    
    ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Method') + ggtitle(title)
  }
}

groupbyCancer(cancer = 1, scoresNormal, title = 'Complete Data')
groupbyCancer(cancer = 2, scoresNormal, title = 'Complete Data')
groupbyCancer(cancer = 3, scoresNormal, title = 'Complete Data')
groupbyCancer(cancer = 4, scoresNormal, title = 'Complete Data')

groupbyCancer(cancer = 1, scoresNormalOrig, orig = TRUE, title = 'Original Data')
groupbyCancer(cancer = 2, scoresNormalOrig, orig = TRUE, title = 'Original Data')
groupbyCancer(cancer = 3, scoresNormalOrig, orig = TRUE, title = 'Original Data')
groupbyCancer(cancer = 4, scoresNormalOrig, orig = TRUE, title = 'Original Data')


###############################################################################################################
# Distribution of acc, nmi, pval, and ci
distComplete <- function(column) {
  par(mfrow = c(2,2))
  hist(scoresNormal[, column], main = paste0('complete data', ' ', column), xlab = 'Normal', col = 'lightblue')
  hist(scoresCombat[, column], main = paste0('complete data', ' ', column), xlab = 'Combat', col = 'lightblue' )
  hist(scoresGenderMale[, column], main = paste0('complete data', ' ', column), xlab = 'Male', col = 'lightblue')
  hist(scoresGenderFemale[, column], main = paste0('complete data', ' ', column), xlab = 'Female', col = 'lightblue')
}
distComplete('acc')
distComplete('nmi')
distComplete('pval')
distComplete('ci')


distOrig <- function(column) {
  par(mfrow =  c(2,2))
  hist(scoresNormalOrig[, column], main = paste0('original data', ' ', column), xlab = 'Normal', col = 'lightblue')
  hist(scoresCombatOrig[, column], main = paste0('original data', ' ', column), xlab = 'Combat', col = 'lightblue')
  hist(scoresGenderOrigMale[, column], main = paste0('original data', ' ', column), xlab = 'Male', col = 'lightblue')
  hist(scoresGenderOrigFemale[, column], main = paste0('original data', ' ', column), xlab = 'Female', col = 'lightblue')
}
distOrig('pval')
distOrig('ci')











