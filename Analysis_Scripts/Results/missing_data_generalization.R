################################################################################################
# This script will take the results from the intersection and union and create a plot that indicates 
# the number of times a method is ranked same in int and union. We want to see which cluster sizes 
# gives bets generalization. 

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
testScores <- read.csv(paste0(results_folder, '/missing_data_int.csv'))
testScoresOrig <- read.csv(paste0(results_folder, '/missing_data_union.csv'))

# remove unneeded columns
testScoresOrig$acc <- NULL
testScoresOrig$nmi <- NULL
testScores$X <- NULL


# get ranking for each method
rankMethods <- function(data, complete) {
  
  if (complete) {
    
    data <- transform(data, 
                      acc_nmi_rank_int = ave(acc_nmi, 
                                         cancer,
                                         cluster,
                                         FUN = function(x) rank(-x, ties.method = "min")),
                      pval_ci_rank_int = ave(pval_ci, 
                                         cancer, 
                                         cluster,
                                         FUN = function(x) rank(-x, ties.method = "min")),
                      total_rank_int = ave(total, 
                                           cancer, 
                                           cluster,
                                           FUN = function(x) rank(-x, ties.method = "min")),
                      pval_rank_int = ave(pval, 
                                      cancer, 
                                      cluster,
                                      FUN = function(x) rank(-x, ties.method = "min")))
    
  } else {
    
    data <- transform(data, 
                      pval_ci_rank_union = ave(total, 
                                             cancer, 
                                             cluster,
                                             FUN = function(x) rank(-x, ties.method = "min")),
                      pval_rank_union = ave(pval, 
                                            cancer,
                                            cluster,
                                            FUN = function(x) rank(-x, ties.method = "min")))
  }
  
  return(data)
}

int_rank <- rankMethods(testScores, complete = TRUE)
union_rank <- rankMethods(testScoresOrig, complete = FALSE)

########################################################################################
# # For total ranking. That is, top ranked for total_rank in intersection and top rank for total_rank in union
# 
# # first assign qunitle ranking to both intersection and union 
# int_rank$quintile_int <- ifelse(int_rank$total_rank > 0 & int_rank$total_rank <= 3, 'first',
#                             ifelse(int_rank$total_rank > 3 & int_rank$total_rank <= 6, 'second',
#                                    ifelse(int_rank$total_rank > 6 & int_rank$total_rank <= 9, 'third',
#                                           ifelse(int_rank$total_rank > 9 & int_rank$total_rank <= 12, 'fourth', 'fifth'))))
# 
# union_rank$quintile_union <- ifelse(union_rank$total_rank > 0 & union_rank$total_rank <= 3, 'first',
#                             ifelse(union_rank$total_rank > 3 & union_rank$total_rank <= 6, 'second',
#                                    ifelse(union_rank$total_rank > 6 & union_rank$total_rank <= 9, 'third',
#                                           ifelse(union_rank$total_rank > 9 & union_rank$total_rank <= 12, 'fourth', 'fifth'))))

rm(testScores, testScoresOrig)
#########################################################################################
## Create data frame that has counts for number of time method is in same quintile in intersection and union 

# first paste cancer and method together 
int_rank$cancer_method_cluster <- paste0(int_rank$cancer,'_', int_rank$method, '_', int_rank$cluster)
union_rank$cancer_method_cluster <- paste0(union_rank$cancer,'_', union_rank$method, '_', union_rank$cluster)

# left join on cancer_method 
dat <- left_join(int_rank, 
                 union_rank, 
                 by = 'cancer_method_cluster')

# remove unnecessary columns 
dat <- dat[, c('cancer.x', 
               'method.x',
               'cluster.x',
               'cancer_method_cluster',
               'pval_rank_int',
               'pval_rank_union',
               'pval_ci_rank_int', 
               'pval_ci_rank_union',
               'acc_nmi_rank_int', 
               'total_rank_int')]

#############################################################################################
## check if method is in within 1 place (up or down) in both int and union 

# for total score of intersection and total union
dat$total <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$total_rank_int[i] == dat$pval_ci_rank_union[i]) {
    dat$total[i] <- TRUE
    
  } else if (abs(dat$total_rank_int[i] - dat$pval_ci_rank_union[i]) == 1) {
    dat$total[i] <- TRUE
    
  } else if (abs(dat$total_rank_int[i] - dat$pval_ci_rank_union[i]) > 1) {
    dat$total[i] <- FALSE
  }
  
}

# for acc_nmi of intersection and total union
dat$total_acc_nmi <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$acc_nmi_rank_int[i] == dat$pval_ci_rank_union[i]) {
    dat$total_acc_nmi[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank_int[i] - dat$pval_ci_rank_union[i]) == 1) {
    dat$total_acc_nmi[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank_int[i] - dat$pval_ci_rank_union[i]) > 1) {
    dat$total_acc_nmi[i] <- FALSE
  }
  
}


# for pval_ci of intersection and total union
dat$total_pval_ci <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$pval_ci_rank_int[i] == dat$pval_ci_rank_union[i]) {
    dat$total_pval_ci[i] <- TRUE
    
  } else if (abs(dat$pval_ci_rank_int[i] - dat$pval_ci_rank_union[i]) == 1) {
    dat$total_pval_ci[i] <- TRUE
    
  } else if (abs(dat$pval_ci_rank_int[i] - dat$pval_ci_rank_union[i]) > 1) {
    dat$total_pval_ci[i] <- FALSE
  }
}

# for individual pval score of intersection and total union
dat$total_pval <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$pval_rank_int[i] == dat$pval_rank_union[i]) {
    dat$total_pval[i] <- TRUE
    
  } else if (abs(dat$pval_rank_int[i] - dat$pval_rank_union[i]) == 1 ||
             abs(dat$pval_rank_int[i] - dat$pval_rank_union[i]) == 2) {
    dat$total_pval[i] <- TRUE
    
  } else if (abs(dat$pval_rank_int[i] - dat$pval_rank_union[i]) > 2) {
    dat$total_pval[i] <- FALSE
  }
  
}


###############################################################################################
# group by cancer 
dat_counts <- dat %>% group_by(cancer.x,method.x) %>% summarise(total = sum(total == TRUE),
                                                       total_acc_nmi = sum(total_acc_nmi == TRUE),
                                                       total_pval_ci = sum(total_pval_ci == TRUE),
                                                       total_pval = sum(total_pval == TRUE))



##############################################################################################
# Find methods that have consistent scores for acc_nmi and pval_ci 

# for total score of intersection and total union
dat$acc_pval <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$acc_nmi_rank_int[i] == dat$pval_ci_rank_int[i]) {
    dat$acc_pval[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank_int[i] - dat$pval_ci_rank_int[i]) == 2 ||
             abs(dat$acc_nmi_rank_int[i] - dat$pval_ci_rank_int[i]) == 1) {
    dat$acc_pval[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank_int[i] - dat$pval_ci_rank_int[i]) > 2) {
    dat$acc_pval[i] <- FALSE
  }
  
}

###############################################################################################
# group by cancer 
dat_counts <- dat %>% group_by(cancer.x, cluster.x) %>% summarise(acc_pval = sum(acc_pval == TRUE))


