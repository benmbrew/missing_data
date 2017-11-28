##### Analyze testScores from LUSC 

######################################################################
# Load libraries
# library(dplyr)
# library(reshape2)

######################################################################
# Initialize folders
home_folder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
project_folder <- paste(home_folder, 'Projects/SNF/NM_2015', sep = '/')
scripts_folder <- paste(project_folder, "Scripts", "06_LUSC", sep = "/")
get_data <- paste(scripts_folder, "analyze_results", sep ="/")

# Load testScores
dat <- read.table(paste(get_data, 'testScores.txt', sep='/'), header = TRUE)
