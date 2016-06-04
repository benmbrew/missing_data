################################################################################################
# This script will create box plots for imputation nrmse

# iniiate folders
################################################################################################
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
################################################################################################
# Initialize folders, 
home_folder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
project_folder <- paste(home_folder, 'Projects/SNF/NM_2015', sep = '/')
scripts_folder <- paste(project_folder, "Scripts", '06_Two_Thousand_Features', sep = "/")
plots_folder <- paste(scripts_folder, "analyze_results/Plots", sep ="/")
results_folder <- paste(project_folder, 'Scripts/06_Results', sep = '/')
source(paste0(results_folder, '/Lib/helpers.R'))

######################################################################
# Initialize fixed variables
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC") 
dataTypes <- c("methyl", "mirna", "mrna")
imputeTypes <- c("KNN", "LLS", "LSA", "Rand")
clusterTypes <- c("Hierarchical", "iCluster", "SNF")
setTypes <- c("Union", "Intersect")
similarityTypes <- c("Self", "Median", "Regres")

######################################################################
# Load the results
impute <- read.table(paste0(results_folder, '/imputation.txt'))

# Load the imputation results
names(impute) <- c("cancer", "impute", "seed", "methyl", "mirna", "mrna",
                "runtime")

#########################################################################
### Create box plot for each method data type 
  
# name variables 
impute$cancer <- cancerTypes[impute$cancer]
impute$impute <- imputeTypes[impute$impute]

# # first create barplot for run time
# ggplot(data = impute, aes(impute, runtime)) + geom_bar(stat = 'identity')

# get rid of runtim
impute$runtime <- NULL
impute <- impute[which(impute$impute != 'LLS'),]

# melt data frame to there is variable that has all data types in grouped form
impute <- melt(impute, id.vars = c('cancer', 'impute', 'seed'), 
             variable.name = 'data_type', 
             value.name = 'nrmse')

# split into cancers 
impute_brca <- impute[impute$cancer == 'BRCA',]
impute_kirc <- impute[impute$cancer == 'KIRC',]
impute_lihc <- impute[impute$cancer == 'LIHC',]
impute_luad <- impute[impute$cancer == 'LUAD',]
impute_lusc <- impute[impute$cancer == 'LUSC',]

# make box plot across all cancers with method on x axis and data type y axis 

# Create a labeller for facet_wrap for box plot 
data_label <- list(
  'methyl' = "Methylation",
  'mirna' = "Mirna",
  'mrna' = "Mrna"
)

data_labeller <- function(variable,value){
  return(data_label[value])
}

# BRCA
ggplot(impute_brca, aes(impute, nrmse)) + 
  geom_boxplot(fill = 'grey50') + facet_wrap(~data_type, labeller = data_labeller) + 
  theme(strip.text.x = element_text(size = 20, colour = "grey30", angle = 0, face = 'bold')) + 
  xlab('Imputation Method') + ylab('nmrse') + ggtitle('BRCA Imputation Scores') +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

# KIRC
ggplot(impute_KIRC, aes(impute, nrmse)) + 
  geom_boxplot(fill = 'grey50') + facet_wrap(~data_type, labeller = data_labeller) + 
  theme(strip.text.x = element_text(size = 20, colour = "grey30", angle = 0, face = 'bold')) + 
  xlab('Imputation Method') + ylab('nmrse') + ggtitle('KIRC Imputation Scores') +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

# LIHC
ggplot(impute_lihc, aes(impute, nrmse)) + 
  geom_boxplot(fill = 'grey50') + facet_wrap(~data_type, labeller = data_labeller) + 
  theme(strip.text.x = element_text(size = 20, colour = "grey30", angle = 0, face = 'bold')) + 
  xlab('Imputation Method') + ylab('nmrse') + ggtitle('LIHC Imputation Scores') +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

# LUAD
ggplot(impute_luad, aes(impute, nrmse)) + 
  geom_boxplot(fill = 'grey50') + facet_wrap(~data_type, labeller = data_labeller) + 
  theme(strip.text.x = element_text(size = 20, colour = "grey30", angle = 0, face = 'bold')) + 
  xlab('Imputation Method') + ylab('nmrse') + ggtitle('LUAD Imputation Scores') +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

# LUSC
ggplot(impute_lusc, aes(impute, nrmse)) + 
  geom_boxplot(fill = 'grey50') + facet_wrap(~data_type, labeller = data_labeller) + 
  theme(strip.text.x = element_text(size = 20, colour = "grey30", angle = 0, face = 'bold')) + 
  xlab('Imputation Method') + ylab('nmrse') + ggtitle('LUSC Imputation Scores') +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 
