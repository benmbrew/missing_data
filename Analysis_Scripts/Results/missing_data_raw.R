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
combat_int <- read.csv(paste0(results_folder, '/missing_data_combat_int.csv'))
combat_union <- read.csv(paste0(results_folder, '/missing_data_combat_union.csv'))

# Load missing_data
int <- read.csv(paste0(results_folder, '/missing_data_complete.csv'))
union <- read.csv(paste0(results_folder, '/missing_data_orig.csv'))

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
int <- int[, c('impute', 'cluster', 'acc', 'nmi', 'pval', 'ci', 'method', 'cancer', 'cluster_size')]
union <- union[, c('impute', 'cluster', 'acc', 'nmi', 'pval', 'ci', 'method', 'cancer', 'cluster_size')]

#######################################################################################################
# subset by cancer and cluster size
# optimal clusters are 3, 2, 4, 4, 4

#intersection
brca_int <- int[int$cancer == 'BRCA' & int$cluster_size == 3,]
kirc_int <- int[int$cancer == 'combat' & int$cluster_size == 3,]
lihc_int <- int[int$cancer == 'LIHC' & int$cluster_size == 4,]
luad_int <- int[int$cancer == 'LUAD' & int$cluster_size == 4,]
lusc_int <- int[int$cancer == 'LUSC' & int$cluster_size == 4,]

# union
brca_union <- union[union$cancer == 'BRCA' & union$cluster_size == 3,]
kirc_union <- union[union$cancer == 'combat' & union$cluster_size == 3,]
lihc_union <- union[union$cancer == 'LIHC' & union$cluster_size == 4,]
luad_union <- union[union$cancer == 'LUAD' & union$cluster_size == 4,]
lusc_union <- union[union$cancer == 'LUSC' & union$cluster_size == 4,]

######################################################################################################
# group by method for each cancer and get mean for intersectio nand union 

groupCancer <- function(data,
                        sig = FALSE,
                        acc = FALSE, 
                        nmi = FALSE, 
                        pval = FALSE, 
                        ci= FALSE, 
                        acc_nmi = FALSE, 
                        pval_ci = FALSE, 
                        ylab,
                        title) {
  if(sig) {
    data <- data[data$pval > 1.3,]
  }
  
  if (acc) {
  group_data <- data %>%
    group_by(method) %>%
    summarise(value = mean(acc, na.rm = T))
  } else if (nmi) {
    group_data <- data %>%
      group_by(method) %>%
      summarise(value = mean(nmi, na.rm = T))
  } else if (pval) {
    group_data <- data %>%
      filter(!is.na(method)) %>%
      group_by(method) %>%
      summarise(value = mean(pval, na.rm = T))
  } else if (ci) {
    group_data <- data %>%
      group_by(method) %>%
      summarise(value = mean(ci, na.rm = T))
  } else if (acc_nmi) {
    group_data <- data %>%
      group_by(method) %>%
      summarise(value = sum(mean(acc, na.rm = T), mean(nmi, na.rm = T)))
  } else if(pval_ci) {
    group_data <- data %>%
      group_by(method) %>%
      summarise(value = sum(mean(pval, na.rm = T), mean(ci, na.rm = T)))
  }
  
  # plot 
  ggplot(data = group_data, aes(reorder(method, -value), value)) +
    geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('cancer') + ylab(ylab) + 
    ggtitle(title) + geom_text(aes(label=round(value , 3)), vjust=1.6, color="white", size=4) + 
    theme(panel.background=element_rect(fill="white"), 
          plot.background=element_rect(fill="white"), 
          panel.border = element_rect(fill = NA, colour = 'grey50'),
          panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
          legend.position="none", #legend.title = element_blank(), 
          legend.background = element_rect(fill="#F0F0F0"),
          plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
          axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
          axis.text.y=element_text(size=11,colour="#535353",face="bold"),
          axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
          axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
          plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)
  
}
# 1.30

# BRCA and intersection
groupCancer(brca_int, acc= TRUE, ylab = 'Acc', title = 'BRCA Intersection')
groupCancer(brca_int, nmi= TRUE, ylab = 'NMI', title = 'BRCA Intersection')
groupCancer(brca_int, pval= TRUE, ylab = 'Pval', title = 'BRCA Intersection')
groupCancer(brca_int, sig = TRUE,  pval= TRUE, ylab = 'Pval', title = 'BRCA Intersection')
groupCancer(brca_int, ci= TRUE, ylab = 'CI', title = 'BRCA Intersection')
groupCancer(brca_int, acc_nmi = TRUE, ylab = 'Acc + NMI', title = 'BRCA Intersection')
groupCancer(brca_int, pval_ci = TRUE, ylab = 'Pval + CI', title = 'BRCA Intersection')

# BRCA and union
groupCancer(brca_union, pval= TRUE, ylab = 'Pval', title = 'BRCA Union')
groupCancer(brca_union, sig = TRUE, pval= TRUE, ylab = 'Pval', title = 'BRCA Union')

groupCancer(brca_union, ci= TRUE, ylab = 'CI', title = 'BRCA Union')
groupCancer(brca_union, pval_ci = TRUE, ylab = 'Pval + CI', title = 'BRCA Union')


# KIRC and intersection
groupCancer(kirc_int, acc= TRUE, ylab = 'Acc', title = 'KIRC Intersection')
groupCancer(kirc_int, nmi= TRUE, ylab = 'NMI', title = 'KIRC Intersection')
groupCancer(kirc_int, pval= TRUE, ylab = 'Pval', title = 'KIRC Intersection')
groupCancer(kirc_int, sig = TRUE, pval= TRUE, ylab = 'Pval', title = 'KIRC Intersection')
groupCancer(kirc_int, ci= TRUE, ylab = 'CI', title = 'KIRC Intersection')
groupCancer(kirc_int, acc_nmi = TRUE, ylab = 'Acc + NMI', title = 'KIRC Intersection')
groupCancer(kirc_int, pval_ci = TRUE, ylab = 'Pval + CI', title = 'KIRC Intersection')

# kirc and union
groupCancer(kirc_union, pval= TRUE, ylab = 'Pval', title = 'KIRC Union')
groupCancer(kirc_union, sig = TRUE,  pval= TRUE, ylab = 'Pval', title = 'KIRC Union')

groupCancer(kirc_union, ci= TRUE, ylab = 'CI', title = 'KIRC Union')
groupCancer(kirc_union, pval_ci = TRUE, ylab = 'Pval + CI', title = 'KIRC Union')

# lihc and intersection
groupCancer(lihc_int, acc= TRUE, ylab = 'Acc', title = 'LIHC Intersection')
groupCancer(lihc_int, nmi= TRUE, ylab = 'NMI', title = 'LIHC Intersection')
groupCancer(lihc_int, pval= TRUE, ylab = 'Pval', title = 'LIHC Intersection')
groupCancer(lihc_int, sig = TRUE, pval= TRUE, ylab = 'Pval', title = 'LIHC Intersection')
groupCancer(lihc_int, ci= TRUE, ylab = 'CI', title = 'LIHC Intersection')
groupCancer(lihc_int, acc_nmi = TRUE, ylab = 'Acc + NMI', title = 'LIHC Intersection')
groupCancer(lihc_int, pval_ci = TRUE, ylab = 'Pval + CI', title = 'LIHC Intersection')

# lihc and union
groupCancer(lihc_union, pval= TRUE, ylab = 'Pval', title = 'LIHC Union')
groupCancer(lihc_union, sig = TRUE, pval= TRUE, ylab = 'Pval', title = 'LIHC Union')
groupCancer(lihc_union, ci= TRUE, ylab = 'CI', title = 'LIHC Union')
groupCancer(lihc_union, pval_ci = TRUE, ylab = 'Pval + CI', title = 'LIHC Union')

# luad and intersection
groupCancer(luad_int, acc= TRUE, ylab = 'Acc', title = 'LUAD Intersection')
groupCancer(luad_int, nmi= TRUE, ylab = 'NMI', title = 'LUAD Intersection')
groupCancer(luad_int, pval= TRUE, ylab = 'Pval', title = 'LUAD Intersection')
groupCancer(luad_int, sig = TRUE, pval= TRUE, ylab = 'Pval', title = 'LUAD Intersection')
groupCancer(luad_int, ci= TRUE, ylab = 'CI', title = 'LUAD Intersection')
groupCancer(luad_int, acc_nmi = TRUE, ylab = 'Acc + NMI', title = 'LUAD Intersection')
groupCancer(luad_int, pval_ci = TRUE, ylab = 'Pval + CI', title = 'LUAD Intersection')

# luad and union
groupCancer(luad_union, pval= TRUE, ylab = 'Pval', title = 'LUAD Union')
groupCancer(luad_union, sig = TRUE, pval= TRUE, ylab = 'Pval', title = 'LUAD Union')
groupCancer(luad_union, ci= TRUE, ylab = 'CI', title = 'LUAD Union')
groupCancer(luad_union, pval_ci = TRUE, ylab = 'Pval + CI', title = 'LUAD Union')

# LUSC and intersection
groupCancer(lusc_int, acc= TRUE, ylab = 'Acc', title = 'LUSC Intersection')
groupCancer(lusc_int, nmi= TRUE, ylab = 'NMI', title = 'LUSC Intersection')
groupCancer(lusc_int, pval= TRUE, ylab = 'Pval', title = 'LUSC Intersection')
groupCancer(lusc_int, ci= TRUE, ylab = 'CI', title = 'LUSC Intersection')
groupCancer(lusc_int, acc_nmi = TRUE, ylab = 'Acc + NMI', title = 'LUSC Intersection')
groupCancer(lusc_int, pval_ci = TRUE, ylab = 'Pval + CI', title = 'LUSC Intersection')

# LUSC and union
groupCancer(lusc_union, pval= TRUE, ylab = 'Pval', title = 'LUSC Union')
groupCancer(lusc_union, ci= TRUE, ylab = 'CI', title = 'LUSC Union')
groupCancer(lusc_union, pval_ci = TRUE, ylab = 'Pval + CI', title = 'LUSC Union')

#######################################################################################
# Look at pvalue for each method and caner 

# seperate data by significant pvalue and examine methods 

pvalByCancer <- function(data, ylab, title) {
  
  data$sig <- ifelse(data$pval > 1.30, TRUE, FALSE)
  
  group <- data %>%
    group_by(method) %>%
    summarise(value = sum(sig == TRUE, na.rm = T))
  
  # plot 
  ggplot(data = group, aes(reorder(method, -value), value)) +
    geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('cancer') + ylab(ylab) + 
    ggtitle(title) + geom_text(aes(label=round(value , 3)), vjust=1.6, color="white", size=4) + 
    theme(panel.background=element_rect(fill="white"), 
          plot.background=element_rect(fill="white"), 
          panel.border = element_rect(fill = NA, colour = 'grey50'),
          panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
          legend.position="none", #legend.title = element_blank(), 
          legend.background = element_rect(fill="#F0F0F0"),
          plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
          axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
          axis.text.y=element_text(size=11,colour="#535353",face="bold"),
          axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
          axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
          plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)
  
}

# for intersection
pvalByCancer(brca_int, ylab = '# of times significant', title = "BRCA Intersection")
pvalByCancer(kirc_int, ylab = '# of times significant', title = "KIRC Intersection")
pvalByCancer(lihc_int, ylab = '# of times significant', title = "LIHC Intersection")
pvalByCancer(luad_int, ylab = '# of times significant', title = "LUAD Intersection")
pvalByCancer(lusc_int, ylab = '# of times significant', title = "LUSC Intersection")

# for Union
pvalByCancer(brca_union, ylab = '# of times significant', title = "BRCA Union")
pvalByCancer(kirc_union, ylab = '# of times significant', title = "KIRC Union")
pvalByCancer(lihc_union, ylab = '# of times significant', title = "LIHC Union")
pvalByCancer(luad_union, ylab = '# of times significant', title = "LUAD Union")
pvalByCancer(lusc_union, ylab = '# of times significant', title = "LUSC Union")
