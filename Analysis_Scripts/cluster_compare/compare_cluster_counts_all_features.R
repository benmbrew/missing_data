####################################################################################################
# This script will look at cluster size for different cluster numbers 
library(ggplot2)
library(reshape2)
library(dplyr)
################################################################################################
# Initialize folders, 
home_folder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
project_folder <- paste(home_folder, 'Projects/SNF/NM_2015', sep = '/')
scripts_folder <- paste(project_folder, "Scripts", '06_Cluster_Sizes_all_features', sep = "/")
plots_folder <- paste(scripts_folder, "analyze_results/Plots", sep ="/")
results_folder <- paste(project_folder, 'Scripts/06_Results', sep = '/')
source(paste0(results_folder, '/Lib/helpers.R'))

####################################################################################################
# Read in all_labels.csv from results folder 
dat <- read.csv(paste0(results_folder, '/all_labels_all_features.csv'))

# get percent 
dat$per_1 <- round((dat$label1/dat$total)*100,2)
dat$per_2 <- round((dat$label2/dat$total)*100,2)
dat$per_3 <- round((dat$label3/dat$total)*100,2)
dat$per_4 <- round((dat$label4/dat$total)*100,2)
dat$per_5 <- round((dat$label5/dat$total)*100,2)

###################################################################################################
# group by cancer and cluster, get counts of label greater than 5

dat_5 <- dat %>%
  group_by(cancer, clusters) %>%
  summarise(total = sum(per_1 < 5, na.rm = T) + sum(per_2 < 5, na.rm = T) + sum(per_3 < 5, na.rm = T) + 
              sum(per_4 < 5, na.rm = T) + sum(per_5 < 5, na.rm = T))

# group by cancer and cluster, get counts of label greater than 10
dat_10 <- dat %>%
  group_by(cancer, clusters) %>%
  summarise(total = sum(per_1 < 10, na.rm = T) + sum(per_2 < 10, na.rm = T) + sum(per_3 < 10, na.rm = T) + 
              sum(per_4 < 10, na.rm = T) + sum(per_5 < 10, na.rm = T))


# group by cancer and cluster, get counts of label greater than 10
dat_20 <- dat %>%
  group_by(cancer, clusters) %>%
  summarise(total = sum(per_1 < 20, na.rm = T) + sum(per_2 < 20, na.rm = T) + sum(per_3 < 20, na.rm = T) + 
              sum(per_4 < 20, na.rm = T) + sum(per_5 < 20, na.rm = T))









