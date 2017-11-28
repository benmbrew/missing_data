######################################################################################
# This script will look at cluster sizes
################################################################################################
# This script will take the results from all folder and make barplots for each strategy and method 
library(ggplot2)
library(reshape2)
library(dplyr)
library(survival)
################################################################################################
# Initialize folders, 
home_folder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
project_folder <- paste(home_folder, 'Projects/SNF/NM_2015', sep = '/')
scripts_folder <- paste(project_folder, "Scripts", '06_Two_Thousand_Features', sep = "/")
plots_folder <- paste(scripts_folder, "analyze_results/Plots", sep ="/")
results_folder <- paste(project_folder, 'Scripts/06_Results', sep = '/')
data_folder <- paste(project_folder,'Data', sep = '/')
source(paste0(results_folder, '/Lib/helpers.R'))


# Load data
all_labels <- read.csv(paste0(results_folder, '/all_lables.csv'))
intersection <- read.csv(paste0(results_folder, '/intersection_clusters.csv'))
union <- read.csv(paste0(results_folder, '/union_clusters.csv'))

##################################################################################################
# clean intersection and union 

# first remove the _clustersiez from cancer column 
intersection$cancer <- gsub('_3', '', intersection$cancer)
intersection$cancer <- gsub('_4', '', intersection$cancer)

union$cancer <- gsub('_3', '', union$cancer)
union$cancer <- gsub('_4', '', union$cancer)

# merge cancer column onto method column
intersection$method <- paste(intersection$cancer, 'com', intersection$method, sep = '_')
union$method <- paste(union$cancer, 'union', union$method, sep = '_')

# add column for union and intersection
intersection$type <- 'intersection'
union$type <- 'union'


# make method column lower case 
intersection$method <- tolower(intersection$method)
union$method <- tolower(union$method)

# replace '.' with "_"
intersection$method <- gsub('.', '_', intersection$method, fixed = TRUE)
union$method <- gsub('.', '_', union$method, fixed = TRUE)

# replace "hierarchical" and "icluster" with "hier" and "iclust"
intersection$method <- gsub('icluster', 'iclust', intersection$method)
union$method <- gsub('icluster', 'iclust', union$method)

intersection$method <- gsub('hierarchical', 'hier', intersection$method)
union$method <- gsub('hierarchical', 'hier', union$method)

intersection$method <- gsub('regres', 'reg', intersection$method)
union$method <- gsub('regres', 'reg', union$method)

intersection$method <- gsub('median', 'med', intersection$method)
union$method <- gsub('median', 'med', union$method)

# paste cluster size at end of method
intersection$method <- paste(intersection$method, intersection$clusters, sep = '_')
union$method <- paste(union$method, union$clusters, sep = '_')

# row bind intersection and union. First need to create NA columns for acc and nmi in union and drop totals
union$acc <- NA
union$nmi <- NA
union$total <- NULL
intersection$total <- NULL
intersection$acc_nmi <- NULL
intersection$pval_ci <- NULL

scores <- rbind(intersection, union)

# remove combat from analysis
scores <- scores[!grepl('combat|Combat', scores$cancer),]


###################################################################################################
# create variables for percent
all_labels$lab1_percent <- round((all_labels$label1/all_labels$sum_of_samples)*100,2)
all_labels$lab2_percent <- round((all_labels$label2/all_labels$sum_of_samples)*100,2)
all_labels$lab3_percent <- round((all_labels$label3/all_labels$sum_of_samples)*100,2)
all_labels$lab4_percent <- round((all_labels$label4/all_labels$sum_of_samples)*100,2)
all_labels$lab5_percent <- round((all_labels$label5/all_labels$sum_of_samples)*100,2)

# add cluster size to end of methods 
all_labels$method <- paste(all_labels$method, all_labels$clusters, sep = '_')

# remove duplicates 
all_labels <- all_labels[!duplicated((all_labels$method)),]

# add in column to indicate if it is complete data 
all_labels$data <- ifelse(nchar(all_labels$method) < 18, 'complete', 'not_complete')

###################################################################################################

# group by method and cluster numbers and get counts for label percent under 1%
percent_1 <- all_labels %>%
  group_by(method, clusters, data) %>%
  summarise(cluster1_under1 = sum(lab1_percent <= 1),
            cluster2_under1 = sum(lab2_percent <= 1),
            cluster3_under1 = sum(lab3_percent <= 1),
            cluster4_under1 = sum(lab4_percent <= 1),
            cluster5_under1 = sum(lab5_percent <= 1))

# left join results from intersection and union onto 
percent_1_scores <- left_join(percent_1, scores)

# fill data= 'notcomplete' with type
notcomplete_index <- percent_1_scores$data == 'not_complete'
percent_1_scores$data[notcomplete_index] <- percent_1_scores$type[notcomplete_index]

# add cancer in column
percent_1_scores$cancer <- ifelse(grepl('brca', percent_1_scores$method), 'BRCA',
                                  ifelse(grepl('kirc', percent_1_scores$method), 'KIRC',
                                         ifelse(grepl('lihc', percent_1_scores$method), 'LIHC',
                                                ifelse(grepl('luad', percent_1_scores$method), 'LUAD',
                                                       ifelse(grepl('lusc', percent_1_scores$method), 
                                                              'LUSC', percent_1_scores$cancer)))))

# group by cancer, get avg acc, nmi, pval, ci and sum of cluster1,2,3,4
group_cluster <- percent_1_scores %>%
  group_by(cancer, data, clusters) %>%
  summarise(total_counts1 = sum(cluster1_under1, na.rm = T),
            total_counts2 = sum(cluster2_under1, na.rm = T),
            total_counts3 = sum(cluster3_under1, na.rm = T),
            total_counts4 = sum(cluster4_under1, na.rm = T),
            total_counts5 = sum(cluster5_under1, na.rm = T))

# add up indicators for all clusters to get a total for under the threshold of x%
group_cluster$total <- group_cluster$total_counts1 + group_cluster$total_counts2 + group_cluster$total_counts3 + 
  + group_cluster$total_counts4 + group_cluster$total_counts5

# write csv 
write.csv(group_cluster, '/home/benbrew/Desktop/one_percent.csv')

#############################################################################################
## group by method and cluster numbers and get counts for label percent under 5%
percent_5 <- all_labels %>%
  group_by(method, clusters, data) %>%
  summarise(cluster1_under5 = sum(lab1_percent <= 5),
            cluster2_under5 = sum(lab2_percent <= 5),
            cluster3_under5 = sum(lab3_percent <= 5),
            cluster4_under5 = sum(lab4_percent <= 5),
            cluster5_under5 = sum(lab5_percent <= 5))

# left join results from intersection and union onto 
percent_5_scores <- left_join(percent_5, scores)

# fill data= 'notcomplete' with type
notcomplete_index <- percent_5_scores$data == 'not_complete'
percent_5_scores$data[notcomplete_index] <- percent_5_scores$type[notcomplete_index]

# add cancer in column
percent_5_scores$cancer <- ifelse(grepl('brca', percent_5_scores$method), 'BRCA',
                                  ifelse(grepl('kirc', percent_5_scores$method), 'KIRC',
                                         ifelse(grepl('lihc', percent_5_scores$method), 'LIHC',
                                                ifelse(grepl('luad', percent_5_scores$method), 'LUAD',
                                                       ifelse(grepl('lusc', percent_5_scores$method), 
                                                              'LUSC', percent_5_scores$cancer)))))

# group by cancer, get avg acc, nmi, pval, ci and sum of cluster5,2,3,4
group_cluster <- percent_5_scores %>%
  group_by(cancer, data, clusters) %>%
  summarise(total_counts1 = sum(cluster1_under5, na.rm = T),
            total_counts2 = sum(cluster2_under5, na.rm = T),
            total_counts3 = sum(cluster3_under5, na.rm = T),
            total_counts4 = sum(cluster4_under5, na.rm = T),
            total_counts5 = sum(cluster5_under5, na.rm = T))

# add up indicators for all clusters to get a total for under the threshold of x%
group_cluster$total <- group_cluster$total_counts1 + group_cluster$total_counts2 + group_cluster$total_counts3 + 
  + group_cluster$total_counts4 + group_cluster$total_counts5

# write csv 
write.csv(group_cluster, '/home/benbrew/Desktop/five_percent.csv')


#####################################################################################################
# group by method and cluster numbers and get counts for label percent under 10%
percent_10 <- all_labels %>%
  group_by(method, clusters, data) %>%
  summarise(cluster1_under10 = sum(lab1_percent <= 10),
            cluster2_under10 = sum(lab2_percent <= 10),
            cluster3_under10 = sum(lab3_percent <= 10),
            cluster4_under10 = sum(lab4_percent <= 10),
            cluster5_under10 = sum(lab5_percent <= 10))

# left join results from intersection and union onto 
percent_10_scores <- left_join(percent_10, scores)

# fill data= 'notcomplete' with type
notcomplete_index <- percent_10_scores$data == 'not_complete'
percent_10_scores$data[notcomplete_index] <- percent_10_scores$type[notcomplete_index]

# add cancer in column
percent_10_scores$cancer <- ifelse(grepl('brca', percent_10_scores$method), 'BRCA',
                                  ifelse(grepl('kirc', percent_10_scores$method), 'KIRC',
                                         ifelse(grepl('lihc', percent_10_scores$method), 'LIHC',
                                                ifelse(grepl('luad', percent_10_scores$method), 'LUAD',
                                                       ifelse(grepl('lusc', percent_10_scores$method), 
                                                              'LUSC', percent_10_scores$cancer)))))

# group by cancer, get avg acc, nmi, pval, ci and sum of cluster10,2,3,4
group_cluster <- percent_10_scores %>%
  group_by(cancer, data, clusters) %>%
  summarise(total_counts1 = sum(cluster1_under10, na.rm = T),
            total_counts2 = sum(cluster2_under10, na.rm = T),
            total_counts3 = sum(cluster3_under10, na.rm = T),
            total_counts4 = sum(cluster4_under10, na.rm = T),
            total_counts5 = sum(cluster5_under10, na.rm = T))

# add up indicators for all clusters to get a total for under the threshold of x%
group_cluster$total <- group_cluster$total_counts1 + group_cluster$total_counts2 + group_cluster$total_counts3 + 
  + group_cluster$total_counts4 + group_cluster$total_counts5

# write csv 
write.csv(group_cluster, '/home/benbrew/Desktop/ten_percent.csv')


#####################################################################################################
# group by method and cluster numbers and get counts for label percent under 15%
percent_15 <- all_labels %>%
  group_by(method, clusters, data) %>%
  summarise(cluster1_under15 = sum(lab1_percent <= 15),
            cluster2_under15 = sum(lab2_percent <= 15),
            cluster3_under15 = sum(lab3_percent <= 15),
            cluster4_under15 = sum(lab4_percent <= 15),
            cluster5_under15 = sum(lab5_percent <= 15))

# left join results from intersection and union onto 
percent_15_scores <- left_join(percent_15, scores)

# fill data= 'notcomplete' with type
notcomplete_index <- percent_15_scores$data == 'not_complete'
percent_15_scores$data[notcomplete_index] <- percent_15_scores$type[notcomplete_index]

# add cancer in column
percent_15_scores$cancer <- ifelse(grepl('brca', percent_15_scores$method), 'BRCA',
                                   ifelse(grepl('kirc', percent_15_scores$method), 'KIRC',
                                          ifelse(grepl('lihc', percent_15_scores$method), 'LIHC',
                                                 ifelse(grepl('luad', percent_15_scores$method), 'LUAD',
                                                        ifelse(grepl('lusc', percent_15_scores$method), 
                                                               'LUSC', percent_15_scores$cancer)))))

# group by cancer, get avg acc, nmi, pval, ci and sum of cluster15,2,3,4
group_cluster <- percent_15_scores %>%
  group_by(cancer, data, clusters) %>%
  summarise(total_counts1 = sum(cluster1_under15, na.rm = T),
            total_counts2 = sum(cluster2_under15, na.rm = T),
            total_counts3 = sum(cluster3_under15, na.rm = T),
            total_counts4 = sum(cluster4_under15, na.rm = T),
            total_counts5 = sum(cluster5_under15, na.rm = T))

# add up indicators for all clusters to get a total for under the threshold of x%
group_cluster$total <- group_cluster$total_counts1 + group_cluster$total_counts2 + group_cluster$total_counts3 + 
  + group_cluster$total_counts4 + group_cluster$total_counts5

# write csv 
write.csv(group_cluster, '/home/benbrew/Desktop/fifteen_percent.csv')


#####################################################################################################
# group by method and cluster numbers and get counts for label percent under 20%
percent_20 <- all_labels %>%
  group_by(method, clusters, data) %>%
  summarise(cluster1_under20 = sum(lab1_percent <= 20),
            cluster2_under20 = sum(lab2_percent <= 20),
            cluster3_under20 = sum(lab3_percent <= 20),
            cluster4_under20 = sum(lab4_percent <= 20),
            cluster5_under20 = sum(lab5_percent <= 20))

# left join results from intersection and union onto 
percent_20_scores <- left_join(percent_20, scores)

# fill data= 'notcomplete' with type
notcomplete_index <- percent_20_scores$data == 'not_complete'
percent_20_scores$data[notcomplete_index] <- percent_20_scores$type[notcomplete_index]

# add cancer in column
percent_20_scores$cancer <- ifelse(grepl('brca', percent_20_scores$method), 'BRCA',
                                   ifelse(grepl('kirc', percent_20_scores$method), 'KIRC',
                                          ifelse(grepl('lihc', percent_20_scores$method), 'LIHC',
                                                 ifelse(grepl('luad', percent_20_scores$method), 'LUAD',
                                                        ifelse(grepl('lusc', percent_20_scores$method), 
                                                               'LUSC', percent_20_scores$cancer)))))

# group by cancer, get avg acc, nmi, pval, ci and sum of cluster20,2,3,4
group_cluster <- percent_20_scores %>%
  group_by(cancer, data, clusters) %>%
  summarise(total_counts1 = sum(cluster1_under20, na.rm = T),
            total_counts2 = sum(cluster2_under20, na.rm = T),
            total_counts3 = sum(cluster3_under20, na.rm = T),
            total_counts4 = sum(cluster4_under20, na.rm = T),
            total_counts5 = sum(cluster5_under20, na.rm = T))

# add up indicators for all clusters to get a total for under the threshold of x%
group_cluster$total <- group_cluster$total_counts1 + group_cluster$total_counts2 + group_cluster$total_counts3 + 
  + group_cluster$total_counts4 + group_cluster$total_counts5


# write csv 
write.csv(group_cluster, '/home/benbrew/Desktop/twenty_percent.csv')

####################################################################################################
# At which threshold do pvalues start to blow up for union data clinical data 

# load clinical data 
brca_clin <- read.table(paste(data_folder,'BRCA_clin.txt', sep = '/'), header = TRUE)
kirc_clin <- read.table(paste(data_folder,'KIRC_clin.txt', sep = '/'), header = TRUE)
lihc_clin <- read.table(paste(data_folder,'LIHC_clin.txt', sep = '/'), header = TRUE)
luad_clin <- read.table(paste(data_folder,'LUAD_clin.txt', sep = '/'), header = TRUE)
lusc_clin <- read.table(paste(data_folder,'LUSC_clin.txt', sep = '/'), header = TRUE)


# define a survival function
survival <- function(clinicalData, 
                     clusterSize, 
                     percent,
                     baseline = FALSE,
                     spread = FALSE) {
  
  if (baseline) {
    remove <- 1/clusterSize
    percent <- remove
    
    label1 <- rep.int(1, percent*nrow(clinicalData))
    label2 <- rep.int(2, remove*nrow(clinicalData))
    label3 <- rep.int(3, remove*nrow(clinicalData))
    label4 <- rep.int(4, remove*nrow(clinicalData))
    label5 <- rep.int(5, remove*nrow(clinicalData))
  } else if (spread) {
    # specify how much to remove after percent is given 
    remove <- abs((1-percent))/(clusterSize-1)
    
    label1 <- rep.int(1, percent*nrow(clinicalData))
    label2 <- rep.int(2, remove*nrow(clinicalData))
    label3 <- rep.int(3, remove*nrow(clinicalData))
    label4 <- rep.int(4, remove*nrow(clinicalData))
    label5 <- rep.int(5, remove*nrow(clinicalData))
    
  } else {
    
    remove <- 1/(clusterSize)
    remove1 <- abs(percent - remove) + remove
    
    label1 <- rep.int(1, percent*nrow(clinicalData))
    label2 <- rep.int(2, remove1*nrow(clinicalData))
    label3 <- rep.int(3, remove*nrow(clinicalData))
    label4 <- rep.int(4, remove*nrow(clinicalData))
    label5 <- rep.int(5, remove*nrow(clinicalData))
  }
  
  if (clusterSize == 5) {
    labels <- c(label1, label2, label3, label4, label5)
  } else if (clusterSize == 4) {
    labels <- c(label1, label2, label3, label4)
    
  } else if (clusterSize == 3) {
    labels <- c(label1, label2, label3)
    
  } else if (clusterSize == 2) {
    labels <- c(label1, label2)
  }
  
  labels <- as.factor(labels)
  
  seed <- 50
  seed_list <- list()
  
  for (i in 1:seed){
    
    set.seed(i)
    labels <- sample(labels)
    
    # subset clinicalData by labels
    clinicalData <- clinicalData[1:length(labels),]
    
    
    survTime <- clinicalData$days_to_death
    deathStatus <- clinicalData$vital_status == "dead"
    missingSurvInd <- is.na(survTime)
    lastFollowup <- clinicalData$days_to_last_followup
    survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
    survObject <- Surv(survTime, deathStatus)
    survDiff <- survdiff(survObject~labels)
    pval <- 1 - pchisq(survDiff$chisq, length(survDiff$n)-1)
    seed_list[i] <- pval
  }
  results <- do.call('rbind', seed_list)
  results <- mean(results)
  return(results)
}

#####################################################################################
# BRCA
# 5 clusters

# baseline is all the clusters spread evenly
survival(brca_clin, clusterSize = 5, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(brca_clin, clusterSize = 5, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(brca_clin, clusterSize = 5, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 5, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 5, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 5, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 5, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 5, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 5, percent = 0.001, baseline = FALSE, spread = FALSE)


# 4 clusters

# baseline is all the clusters spread evenly
survival(brca_clin, clusterSize = 4, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(brca_clin, clusterSize = 4, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(brca_clin, clusterSize = 4, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 4, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 4, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 4, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 4, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 4, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 4, percent = 0.001, baseline = FALSE, spread = FALSE)


# 3 clusters

# baseline is all the clusters spread evenly
survival(brca_clin, clusterSize = 3, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(brca_clin, clusterSize = 3, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(brca_clin, clusterSize = 3, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 3, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 3, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 3, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 3, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 3, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 3, percent = 0.001, baseline = FALSE, spread = FALSE)


# 2 clusters

# baseline is all the clusters spread evenly
survival(brca_clin, clusterSize = 2, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(brca_clin, clusterSize = 2, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(brca_clin, clusterSize = 2, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 2, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 2, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 2, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 2, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(brca_clin, clusterSize = 2, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(brca_clin, clusterSize = 2, percent = 0.001, baseline = FALSE, spread = FALSE)

#####################################################################################
# kirc
# 5 clusters

# baseline is all the clusters spread evenly
survival(kirc_clin, clusterSize = 5, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(kirc_clin, clusterSize = 5, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(kirc_clin, clusterSize = 5, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 5, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 5, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 5, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 5, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 5, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 5, percent = 0.001, baseline = FALSE, spread = FALSE)


# 4 clusters

# baseline is all the clusters spread evenly
survival(kirc_clin, clusterSize = 4, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(kirc_clin, clusterSize = 4, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(kirc_clin, clusterSize = 4, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 4, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 4, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 4, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 4, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 4, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 4, percent = 0.001, baseline = FALSE, spread = FALSE)


# 3 clusters

# baseline is all the clusters spread evenly
survival(kirc_clin, clusterSize = 3, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(kirc_clin, clusterSize = 3, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(kirc_clin, clusterSize = 3, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 3, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 3, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 3, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 3, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 3, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 3, percent = 0.001, baseline = FALSE, spread = FALSE)


# 2 clusters

# baseline is all the clusters spread evenly
survival(kirc_clin, clusterSize = 2, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(kirc_clin, clusterSize = 2, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(kirc_clin, clusterSize = 2, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 2, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 2, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 2, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 2, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(kirc_clin, clusterSize = 2, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(kirc_clin, clusterSize = 2, percent = 0.001, baseline = FALSE, spread = FALSE)

#####################################################################################
# lihc
# 5 clusters

# baseline is all the clusters spread evenly
survival(lihc_clin, clusterSize = 5, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(lihc_clin, clusterSize = 5, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(lihc_clin, clusterSize = 5, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 5, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 5, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 5, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 5, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 5, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 5, percent = 0.001, baseline = FALSE, spread = FALSE)


# 4 clusters

# baseline is all the clusters spread evenly
survival(lihc_clin, clusterSize = 4, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(lihc_clin, clusterSize = 4, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(lihc_clin, clusterSize = 4, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 4, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 4, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 4, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 4, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 4, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 4, percent = 0.001, baseline = FALSE, spread = FALSE)


# 3 clusters

# baseline is all the clusters spread evenly
survival(lihc_clin, clusterSize = 3, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(lihc_clin, clusterSize = 3, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(lihc_clin, clusterSize = 3, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 3, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 3, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 3, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 3, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 3, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 3, percent = 0.001, baseline = FALSE, spread = FALSE)


# 2 clusters

# baseline is all the clusters spread evenly
survival(lihc_clin, clusterSize = 2, baseline = TRUE)

# this has random removal of one cluster, rest spread evenly
survival(lihc_clin, clusterSize = 2, percent = 0.15, baseline = FALSE, spread = TRUE)
# random smaller cluster, one cluster takes all, rest spread evenly
survival(lihc_clin, clusterSize = 2, percent = 0.15, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 2, percent = 0.1, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 2, percent = 0.1, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 2, percent = 0.05, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 2, percent = 0.05, baseline = FALSE, spread = FALSE)

survival(lihc_clin, clusterSize = 2, percent = 0.001, baseline = FALSE, spread = TRUE)
survival(lihc_clin, clusterSize = 2, percent = 0.001, baseline = FALSE, spread = FALSE)













































