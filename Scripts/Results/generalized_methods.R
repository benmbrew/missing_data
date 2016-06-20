################################################################################################
# This script will take the results from the intersection and union and create a plot that indicates 
# the number of times a method is ranked in the same quintile for intersection and
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
testScores <- read.csv(paste0(results_folder, '/testScores.csv'))
testScoresOrig <- read.csv(paste0(results_folder, '/testScoresOrig.csv'))

# remove unneeded columns
testScoresOrig$acc <- NULL
testScoresOrig$nmi <- NULL
testScoresOrig$clusters <- NULL
testScores$union_rank <- NULL
testScores$cluster <- NULL

# add acc and nmi, pval and ci
testScores$acc_nmi <- testScores$acc + testScores$nmi
testScores$pval_ci <- testScores$pval + testScores$ci
testScores$total <- testScores$acc_nmi + testScores$pval_ci

testScoresOrig$total <- testScoresOrig$pval + testScoresOrig$ci

# subset data to just look at SNF
testScores <- testScores[grepl('SNF', testScores$method),]
testScoresOrig <- testScoresOrig[grepl('SNF', testScoresOrig$method),]


# get ranking for each method
rankMethods <- function(data, complete) {
  
  if (complete) {
    
    data <- transform(data, 
                      acc_nmi_rank = ave(acc_nmi, 
                                         cancer, 
                                         FUN = function(x) rank(-x, ties.method = "min")),
                      pval_ci_rank = ave(pval_ci, 
                                         cancer, 
                                         FUN = function(x) rank(-x, ties.method = "min")),
                      total_rank_int = ave(total, 
                                       cancer, 
                                       FUN = function(x) rank(-x, ties.method = "min")),
                      pval_rank = ave(pval, 
                                       cancer, 
                                       FUN = function(x) rank(-x, ties.method = "min")))
    
  } else {
    
    data <- transform(data, 
                      total_rank_union = ave(total, 
                                       cancer, 
                                       FUN = function(x) rank(-x, ties.method = "min")),
                      pval_rank_union = ave(pval, 
                                      cancer, 
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
int_rank$cancer_method <- paste0(int_rank$cancer,'_', int_rank$method)
union_rank$cancer_method <- paste0(union_rank$cancer,'_', union_rank$method)

# left join on cancer_method 
dat <- left_join(int_rank, 
                 union_rank, 
                 by = 'cancer_method')

# remove unnecessary columns 
dat <- dat[, c('cancer.x', 
               'method.x',
               'cancer_method',
               'pval_rank',
               'pval_rank_union',
               'pval_ci_rank', 
               'acc_nmi_rank', 
               'total_rank_int', 
               'total_rank_union')]

#############################################################################################
## check if method is in within 1 place (up or down) in both int and union 

# for total score of intersection and total union
dat$total <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$total_rank_int[i] == dat$total_rank_union[i]) {
    dat$total[i] <- TRUE
    
  } else if (abs(dat$total_rank_int[i] - dat$total_rank_union[i]) == 1 ||
             abs(dat$total_rank_int[i] - dat$total_rank_union[i]) == 2 ||
             abs(dat$total_rank_int[i] - dat$total_rank_union[i]) == 3) {
    dat$total[i] <- TRUE
    
  } else if (abs(dat$total_rank_int[i] - dat$total_rank_union[i]) > 3) {
    dat$total[i] <- FALSE
  }
  
}

# for acc_nmi of intersection and total union
dat$total_acc_nmi <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$acc_nmi_rank[i] == dat$total_rank_union[i]) {
    dat$total_acc_nmi[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank[i] - dat$total_rank_union[i]) == 1 ||
             abs(dat$acc_nmi_rank[i] - dat$total_rank_union[i]) == 2 ||
             abs(dat$acc_nmi_rank[i] - dat$total_rank_union[i]) == 3) {
    dat$total_acc_nmi[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank[i] - dat$total_rank_union[i]) > 3) {
    dat$total_acc_nmi[i] <- FALSE
  }
  
}


# for pval_ci of intersection and total union
dat$total_pval_ci <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$pval_ci_rank[i] == dat$total_rank_union[i]) {
    dat$total_pval_ci[i] <- TRUE
    
  } else if (abs(dat$pval_ci_rank[i] - dat$total_rank_union[i]) == 1 ||
             abs(dat$pval_ci_rank[i] - dat$total_rank_union[i]) == 2 ||
             abs(dat$pval_ci_rank[i] - dat$total_rank_union[i]) == 3) {
    dat$total_pval_ci[i] <- TRUE
    
  } else if (abs(dat$pval_ci_rank[i] - dat$total_rank_union[i]) > 3) {
    dat$total_pval_ci[i] <- FALSE
  }
  
}

# for individual pval score of intersection and total union
dat$total_pval <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$pval_rank[i] == dat$pval_rank_union[i]) {
    dat$total_pval[i] <- TRUE
    
  } else if (abs(dat$pval_rank[i] - dat$pval_rank_union[i]) == 1 ||
             abs(dat$pval_rank[i] - dat$pval_rank_union[i]) == 2 ||
             abs(dat$pval_rank[i] - dat$pval_rank_union[i]) == 3) {
    dat$total_pval[i] <- TRUE
    
  } else if (abs(dat$pval_rank[i] - dat$pval_rank_union[i]) > 3) {
    dat$total_pval[i] <- FALSE
  }
  
}

###############################################################################################

# grab top five for each cancer 
temp <- arrange(dat, cancer.x, total_rank_int) %>%
  group_by(cancer.x) %>%
  summarise(total = sum(total[1:15] == TRUE)) 

# group by cancer 
dat_count <- dat %>%
  group_by(cancer.x) %>% 
  summarise(total = sum(total == TRUE),
            total_acc_nmi = sum(total_acc_nmi == TRUE),
            total_pval_ci = sum(total_pval_ci == TRUE),
            total_pval = sum(total_pval == TRUE))

# group by method 
dat_count <- dat %>%
  group_by(method.x) %>% 
  summarise(total = sum(total == TRUE),
            total_acc_nmi = sum(total_acc_nmi == TRUE),
            total_pval_ci = sum(total_pval_ci == TRUE),
            total_pval = sum(total_pval == TRUE))


# group by cluster method
temp <- strsplit(as.character(dat$method.x), '.', fixed = TRUE)
temp1 <- do.call('rbind', temp)
dat$cluster <- as.factor(temp1[,1])

# group by method 
dat_count <- dat %>%
  group_by(cluster) %>% 
  summarise(total = sum(total == TRUE),
            total_acc_nmi = sum(total_acc_nmi == TRUE),
            total_pval_ci = sum(total_pval_ci == TRUE),
            total_pval = sum(total_pval == TRUE))


##############################################################################################
# Find methods that have consistent scores for acc_nmi and pval_ci 

# for total score of intersection and total union
dat$acc_pval <- NA
for (i in 1:nrow(dat)) {
  
  if (dat$acc_nmi_rank[i] == dat$pval_ci_rank[i]) {
    dat$acc_pval[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank[i] - dat$pval_ci_rank[i]) == 2 ||
             abs(dat$acc_nmi_rank[i] - dat$pval_ci_rank[i]) == 1) {
    dat$acc_pval[i] <- TRUE
    
  } else if (abs(dat$acc_nmi_rank[i] - dat$pval_ci_rank[i]) > 2) {
    dat$acc_pval[i] <- FALSE
  }
  
}

###############################################################################################
# group by cancer 
dat_counts <- dat %>% 
  group_by(cancer.x) %>% 
  summarise(acc_pval = sum(acc_pval == TRUE))

# group by method
dat_counts <- dat %>% 
  group_by(method.x) %>% 
  summarise(acc_pval = sum(acc_pval == TRUE))

# group by cluster_method
dat_counts <- dat %>% 
  group_by(cluster) %>% 
  summarise(acc_pval = sum(acc_pval == TRUE))

