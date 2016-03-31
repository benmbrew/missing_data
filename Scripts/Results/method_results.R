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
scoresNormalOrigClust <- read.csv(paste0(results_folder, '/scoresTwoThousandOrigClust.csv'))


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
scoresNormalOrigInt$method <- interaction(scoresNormalOrigInt$cluster,
                                       scoresNormalOrigInt$impute, drop = TRUE)
scoresNormalOrigClust$method <- interaction(scoresNormalOrigClust$cluster,
                                       scoresNormalOrigClust$impute, drop = TRUE)

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
scoresNormalOrigClust$acc <- NULL
scoresNormalOrigClust$nmi <- NULL
scoresCombatOrig$acc <- NULL
scoresCombatOrig$nmi <- NULL
scoresGenderOrig$acc <- NULL
scoresGenderOrig$nmi <- NULL

####################################################################################################
# Group by method and get mean for each evaluation method (acc, nmi, pval, ci)

groupbyMethod <- function(data, orig = FALSE, title) {
  
   if (orig) {
      temp <- data %>%
        group_by(method) %>%
        summarise(meanPval = mean(pval, na.rm = T),
                  meanCi = mean(ci, na.rm = T))
      
      temp_melt <- melt(temp, id.vars = c('method'))
      ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
         xlab('Method') + ggtitle(title) + theme_538_bar
      
  } else {
    
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


groupbyMethod(scoresNormal, title = 'Complete Data')
groupbyMethod(scoresNormalOrig, orig = TRUE, title = 'Original Data')
groupbyMethod(scoresNormalOrigClust, orig = TRUE, title = 'Original Data')


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

groupbyCancer(cancer = 1, scoresNormal, title = 'BRCA Complete Data')
groupbyCancer(cancer = 2, scoresNormal, title = 'KIRC Complete Data')
groupbyCancer(cancer = 3, scoresNormal, title = 'LIHC Complete Data')
groupbyCancer(cancer = 4, scoresNormal, title = 'LUAD Complete Data')

groupbyCancer(cancer = 1, scoresNormalOrig, orig = TRUE, title = 'BRCA Original Data')
groupbyCancer(cancer = 2, scoresNormalOrig, orig = TRUE, title = 'KIRC Original Data')
groupbyCancer(cancer = 3, scoresNormalOrig, orig = TRUE, title = 'LIHC Original Data')
groupbyCancer(cancer = 4, scoresNormalOrig, orig = TRUE, title = 'LUAD Original Data')


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

groupbyPval(scoresNormalOrigClust, title = 'Union Data BRCA Clust', cancer = 1)
groupbyPval(scoresNormalOrigClust, title = 'Union Data KIRC Clust', cancer = 2)
groupbyPval(scoresNormalOrigClust, title = 'Union Data LIHC Clust', cancer = 3)
groupbyPval(scoresNormalOrigClust, title = 'Union Data LUAD Clust', cancer = 4)


groupbyPval(scoresNormalOrig, title = 'Original Data')
groupbyPval(scoresNormalOrigClust, title = 'Original Data')


groupbyPval(scoresCombat, title = 'KIRC Complete with Combat')
groupbyPval(scoresCombatOrig, title = 'KIRC Original with Combat')
groupbyPval(scoresCombat, scoresNormal, data_indicator = 'Combat', title = 'KIRC Complete with Combat and Normal')
groupbyPval(scoresCombatOrig, scoresNormalOrig, data_indicator = 'Combat', title = 'KIRC Original with Combat and Normal')

groupbyPval(scoresGenderMale, title = 'Complete Data Male')
groupbyPval(scoresGenderOrigMale, title = 'Original Data Male')

groupbyPval(scoresGenderFemale, title = 'Complete Data with Female')
groupbyPval(scoresGenderOrigFemale,title = 'Original Data with Female')

################################################################################################################
# Look at intersection
temp <- scoresNormalOrigInt %>%
  group_by(method) %>%
  summarise(meanAcc = mean(intAcc, na.rm = T),
            meanNmi = mean(intNmi, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
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

######################################################################################################################
# different measure of pvalue and ci

########################################################################################
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

##########################################################################################################
# Original data

scoresNormalOrigPval <- read.csv(paste0(results_folder, '/scoresTwoThousandPval.csv'))

scoresNormalOrigPval$method <- interaction(scoresNormalOrigPval$cluster,
                                       scoresNormalOrigPval$impute, drop = TRUE)

# across cancer pvalcox
temp <- scoresNormalOrigPval %>%
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
temp <- scoresNormalOrigPval %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanCon_index_ci = mean(con_index_ci, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index_ci') + theme_538_bar


# across cancer bias_corrected_c_i
temp <- scoresNormalOrigPval %>%
  group_by(method) %>%
  summarise(meanAcc = mean(acc, na.rm = T),
            meanNmi = mean(nmi, na.rm = T),
            meanPval = mean(pval, na.rm = T),
            meanBias_corrected_c_index = mean(bias_corrected_c_index, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('con_index') + theme_538_bar

# across cancer coef
temp <- scoresNormalOrigPval %>%
  group_by(method) %>%
  summarise(meancoef = mean(coefcox, na.rm = T),
            meanse.std = mean(se.std, na.rm = T))

temp_melt <- melt(temp, id.vars = c('method'))

ggplot(data = temp_melt, aes(reorder(method, -value), value, fill = variable)) + geom_bar(stat = 'identity') +
  xlab('Method') + ggtitle('coef') + theme_538_bar


########################################################################## BRCA 
brca <- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 1,]
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

kirc <- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 2,]
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
lihc <- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 3,]
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

luad<- scoresNormalOrigPval[scoresNormalOrigPval$cancer == 4,]
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
