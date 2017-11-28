
######################################################## 
# Script for exaining the test_scores for 2000 features

########################################################
# Load data 
library(tidyr)
library(readr)
library(ggplot2)
library(dplyr)
library(reshape)
library(splitstackshape)
library(tableplot)

data <- '/home/benbrew/Documents/NM_2015_ben/Scripts/06_two_thousand_features/analyze_results/'
setwd(data)

# load data from 2000 features
test_scores <- read.table('testScores.txt', header = TRUE)

# which method has the most points 
test_scores$acc_nmi = test_scores$acc + test_scores$nmi
test_scores$surv = test_scores$pval + test_scores$ci

#extract max for each cancer type for survival and accuracy measurements
maxMethod <- function(data){
  
  data$cancer_num <- as.numeric(data$cancer)
  measure <- list()
  surv <- list()
  cancer_types <- unique(data$cancer)
  
  for(i in data$cancer_num){
    temp_cancer <- data[data$cancer_num == i,]
    measure[[i]] <- temp_cancer[temp_cancer$acc_nmi == max(temp_cancer$acc_nmi),]
    surv[[i]] <- temp_cancer[temp_cancer$surv == max(temp_cancer$surv),]
  }
  
  
  measure <- do.call('rbind', measure)
  measure$eval_indicator <- 'acc_nmi'
  surv <- do.call('rbind', surv)
  surv$eval_indicator <- 'surv'
  score_result <- rbind(measure, surv)
  
  
  return(score_result)
}

top_scores <- maxMethod(test_scores)

# split into two data sets by eval_indicator
top_scores_acc <- top_scores[top_scores$eval_indicator == 'acc_nmi',]
top_scores_surv <- top_scores[top_scores$eval_indicator == 'surv',]

######################################################################
# Examine original data rank results 

# load data from 2000 features
test_scores_original <- read.table('testScoresOriginal.txt', header = TRUE)

# which method has the most points 
test_scores_original$acc_nmi = test_scores_original$acc + test_scores_original$nmi
test_scores_original$surv = test_scores_original$pval + test_scores_original$ci

#extract max for each cancer type for survival and accuracy measurements
maxMethod <- function(data){
  
  data$cancer_num <- as.numeric(data$cancer)
  measure <- list()
  surv <- list()
  cancer_types <- unique(data$cancer)
  
  for(i in data$cancer_num){
    temp_cancer <- data[data$cancer_num == i,]
    measure[[i]] <- temp_cancer[temp_cancer$acc_nmi == max(temp_cancer$acc_nmi),]
    surv[[i]] <- temp_cancer[temp_cancer$surv == max(temp_cancer$surv),]
  }
  
  
  measure <- do.call('rbind', measure)
  measure$eval_indicator <- 'acc_nmi'
  surv <- do.call('rbind', surv)
  surv$eval_indicator <- 'surv'
  score_result <- rbind(measure, surv)
  
  
  return(score_result)
}

top_scores_original <- maxMethod(test_scores_original)

# split into two data sets by eval_indicator
top_scores_original_acc <- top_scores_original[top_scores_original$eval_indicator == 'acc_nmi',]
top_scores_original_surv <- top_scores_original[top_scores_original$eval_indicator == 'surv',]



