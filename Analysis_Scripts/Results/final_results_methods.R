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
int <- read.csv(paste0(results_folder, '/missing_data_complete.csv'))
union <- read.csv(paste0(results_folder, '/missing_data_orig.csv'))
combat_int <- read.csv(paste0(results_folder, '/missing_data_combat_int.csv'))
combat_union <- read.csv(paste0(results_folder, '/missing_data_combat_union.csv'))

cancerTypes <- c('BRCA', 'KIRC', 'LIHC', 'LUAD', 'LUSC')

###################################################################################################
# Get method types
int$method <- interaction(int$cluster,
                            int$impute, drop = TRUE)

union$method <- interaction(union$cluster,
                          union$impute, drop = TRUE)

combat_int$method <- interaction(combat_int$cluster,
                          combat_int$impute, drop = TRUE)

combat_union$method <- interaction(combat_union$cluster,
                            combat_union$impute, drop = TRUE)

##########################################################################################################
# combine combat with regular data 

# put name of cancer to replace numbers
int$cancer <- cancerTypes[int$cancer]
union$cancer <- cancerTypes[union$cancer]

# create new column
combat_int$cancer <- 'KIRC_Combat'
combat_union$cancer <- 'KIRC_Combat'

# rbind 
int <- rbind(int, combat_int)
union <- rbind(union, combat_union)

# remove KIRC
int <- int[int$cancer != 'KIRC',]
union <- union[union$cancer != 'KIRC',]

####################################################################################################################
# graphs with error bars
groupError <- function(data, 
                       normal = FALSE,
                       by_cancer = FALSE,
                       cancer = FALSE,
                       orig = FALSE, 
                       pval = FALSE,
                       all = FALSE,
                       acc_nmi = FALSE,
                       acc = FALSE,
                       nmi = FALSE, 
                       title) {
  if (by_cancer) {
  data <- data[data$cancer == cancer,]
  } else {
    data <- data
  }
  
  if (orig && pval) {
    temp <- data %>%
      group_by(method) %>%
      summarise(meanPval = mean(pval, na.rm = T))
    
    ggplot(data = temp, aes(reorder(method, -meanPval), meanPval)) + geom_bar(stat = 'identity') +
      xlab('Method') + ylab('Mean -log10 Pval') + ggtitle(title) + theme_538_bar
    
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
      xlab('Method') + ylab('Mean -log10 Pval') + ggtitle(title) + 
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
      xlab('Method') + ylab('Mean NMI') + ggtitle(title) + 
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

groupError(int, normal = TRUE, acc_nmi = TRUE, title = 'Intersection All Cancers')
groupError(int, normal = TRUE, acc = TRUE, title = 'Intersection All Cancers')
groupError(int, normal = TRUE, pval = TRUE, title = 'Intersection All Cancers')
groupError(int, normal = TRUE, nmi = TRUE, title = 'Intersection All Cancers')

groupError(union, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers')

##################################################################################
# by cancer
groupError(int, by_cancer = TRUE, cancer = 'BRCA', normal = TRUE, acc_nmi = TRUE, title = 'Intersection BRCA')
groupError(union, by_cancer = TRUE, cancer = 'BRCA', orig = TRUE, pval = TRUE, title = 'Pval Union BRCA')

groupError(int, by_cancer = TRUE, cancer = 'KIRC_Combat', normal = TRUE, acc = TRUE, title = 'Intersection KIRC')
groupError(union, by_cancer = TRUE, cancer = 'KIRC_Combat', orig = TRUE, pval = TRUE, title = 'Pval Union KIRC')

groupError(int, by_cancer = TRUE, cancer = 'LIHC', normal = TRUE, pval = TRUE, title = 'Intersection LIHC')
groupError(union, by_cancer = TRUE, cancer = 'LIHC', orig = TRUE, pval = TRUE, title = 'Pval Union LIHC')

groupError(int, by_cancer = TRUE, cancer = 'LUAD', normal = TRUE, nmi = TRUE, title = 'Intersection LUAD')
groupError(union, by_cancer = TRUE, cancer = 'LUAD', orig = TRUE, pval = TRUE, title = 'Pval Union LUAD')

groupError(int, by_cancer = TRUE, cancer = 'LUSC', normal = TRUE, nmi = TRUE, title = 'Intersection LUSC')
groupError(union, by_cancer = TRUE, cancer = 'LUAD', orig = TRUE, pval = TRUE, title = 'Pval Union LUSC')


##################################################################################
# remove pvalues in both in and union and combat that are not significant (less than 1.30103)
int_sub <- int[int$pval > 1.30103,]
union_sub <- union[union$pval > 1.30103,]


groupError(int_sub, normal = TRUE, acc_nmi = TRUE, title = 'Intersection All Cancers')
groupError(int_sub, normal = TRUE, acc = TRUE, title = 'Intersection All Cancers')
groupError(int_sub, normal = TRUE, pval = TRUE, title = 'Intersection All Cancers')
groupError(int_sub, normal = TRUE, nmi = TRUE, title = 'Intersection All Cancers')

groupError(union_sub, orig = TRUE, pval = TRUE, title = 'Pval Union All Cancers')


