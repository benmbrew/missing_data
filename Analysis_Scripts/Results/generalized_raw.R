# This script will look at raw scores from final clusters in the intersection and union (with combat) and compate.library(ggplot2)
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

# Set fixed variables
cancerTypes <- c('BRCA', 'KIRC', 'LIHC', 'LUAD', 'LUSC')

# Load combat data 
combat_int <- read.csv(paste0(results_folder, '/scoresCombatDup3.csv'))
combat_union <- read.csv(paste0(results_folder, '/scoresCombatOrigDup3.csv'))
                       
# Load missing_data
int <- read.csv(paste0(results_folder, '/scoresTwoThousandDupClust3.csv'))
union <- read.csv(paste0(results_folder, '/scoresOrigTwoThousandDupClust3.csv'))

# add in cancer type for each data frame
combat_int$cancer <- 'combat'
combat_union$cancer <- 'combat'

int$cancer <- cancerTypes[int$cancer]
union$cancer <- cancerTypes[union$cancer]


# add in method for each data frame 
combat_int$method <- interaction(combat_int$cluster, 
                                 combat_int$impute, 
                                 drop = T)

combat_union$method <- interaction(combat_union$cluster, 
                                   combat_union$impute, 
                                   drop = T)

int$method <- interaction(int$cluster,
                          int$impute,
                          drop = T)

union$method <- interaction(union$cluster,
                          union$impute,
                          drop = T)


# make columns the same for combat_int and int, and combat_union  and union. 
int <- rbind(combat_int, int)
union <- rbind(combat_union,union)
int$X <- NULL
union$X <- NULL

# rename acc, nmi, pval, ci
int <- int[, c('impute', 'cluster', 'acc', 'nmi', 'pval', 'ci', 'method', 'cancer')]
union <- union[, c('impute', 'cluster', 'acc', 'nmi', 'pval', 'ci', 'method', 'cancer')]

names(int) <- c('impute', 'cluster', 'acc_int', 'nmi_int', 'pval_int', 'ci_int', 'method', 'cancer')
names(union) <- c('impute', 'cluster', 'acc_union', 'nmi_union', 'pval_union', 'ci_union', 'method', 'cancer')

# get mu for each cancer and method 
union$cancer_method <- paste0(union$cancer, '_', union$method)
int$cancer_method <- paste0(int$cancer, '_', int$method)

# left_join 
dat <- left_join(int, union, by = 'cancer_method')

# subset dat 
dat <- dat[grepl('SNF', dat$method.x),]

###########################################################################################3
# for pval
# group by cancer and get avg pval for int and union
cancer <- dat %>%
  group_by(cancer.x) %>%
  summarise(pval_int = mean(pval_int, na.rm = T),
            pval_union = mean(pval_union, na.rm = T))

cancer$diff <- abs(cancer$pval_int - cancer$pval_union)

# group by method and get avg pval for int and union
method <- dat %>%
  group_by(method.x) %>%
  summarise(pval_int = mean(pval_int, na.rm = T),
            pval_union = mean(pval_union, na.rm = T))

method$diff <- abs(method$pval_int - method$pval_union)

# group by cancer and method get avg pval for int and union
cancer_method <- dat %>%
  group_by(cancer.x, method.x) %>%
  summarise(pval_int = mean(pval_int, na.rm = T),
            pval_union = mean(pval_union, na.rm = T))

cancer_method$diff <- abs(cancer_method$pval_int - cancer_method$pval_union)


##############################################################################
# for ci

# group by cancer and get avg ci for int and union
cancer <- dat %>%
  group_by(cancer.x) %>%
  summarise(ci_int = mean(ci_int, na.rm = T),
            ci_union = mean(ci_union, na.rm = T))

cancer$diff <- abs(cancer$ci_int - cancer$ci_union)

# group by method and get avg ci for int and union
method <- dat %>%
  group_by(method.x) %>%
  summarise(ci_int = mean(ci_int, na.rm = T),
            ci_union = mean(ci_union, na.rm = T))

method$diff <- abs(method$ci_int - method$ci_union)

# group by cancer and method get avg ci for int and union
cancer_method <- dat %>%
  group_by(cancer.x, method.x) %>%
  summarise(ci_int = mean(ci_int, na.rm = T),
            ci_union = mean(ci_union, na.rm = T))

cancer_method$diff <- abs(cancer_method$ci_int - cancer_method$ci_union)

