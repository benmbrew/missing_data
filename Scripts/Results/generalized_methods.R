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

# add acc and nmi, pval and ci
testScores$acc_nmi <- testScores$acc + testScores$nmi
testScores$pval_ci <- testScores$pval + testScores$ci
testScores$total <- testScores$acc_nmi + testScores$pval_ci

testScoresOrig$total <- testScoresOrig$pval + testScoresOrig$ci


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
                      total_rank = ave(total, 
                                       cancer, 
                                       FUN = function(x) rank(-x, ties.method = "min")))
    
  } else {
    
    data <- transform(data, 
                      total_rank = ave(total, 
                                       cancer, 
                                       FUN = function(x) rank(-x, ties.method = "min")))
  }
  
  return(data)
}

int_rank <- rankMethods(testScores, complete = TRUE)
union_rank <- rankMethods(testScoresOrig, complete = FALSE)

########################################################################################
# For total ranking. That is, top ranked for total_rank in intersection and top rank for total_rank in union

# first assign qunitle ranking to both intersection and union 
int_rank$quintile <- ifelse(int_rank$total_rank > 0 & int_rank$total_rank <= 3, 'first',
                            ifelse(int_rank$total_rank > 3 & int_rank$total_rank <= 6, 'second',
                                   ifelse(int_rank$total_rank > 6 & int_rank$total_rank <= 9, 'third',
                                          ifelse(int_rank$total_rank > 9 & int_rank$total_rank <= 12, 'fourth', 'fifth'))))

union_rank$quunionile <- ifelse(union_rank$total_rank > 0 & union_rank$total_rank <= 3, 'first',
                            ifelse(union_rank$total_rank > 3 & union_rank$total_rank <= 6, 'second',
                                   ifelse(union_rank$total_rank > 6 & union_rank$total_rank <= 9, 'third',
                                          ifelse(union_rank$total_rank > 9 & union_rank$total_rank <= 12, 'fourth', 'fifth'))))
