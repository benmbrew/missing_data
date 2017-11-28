#########################################################################################################
# This script will compare labels from the intersection with those from the union with all features

# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)
library(dplyr)
library(sva)

#########################################################################################################
# initiate folders
homeFolder <- "/home/benbrew/hpf/largeprojects/agoldenb/ben"
projectFolder <- paste(homeFolder, "Projects/SNF/NM_2015", sep="/")
dataFolder <- paste(projectFolder, 'Data', sep = '/')
testFolder <- paste(projectFolder, "Scripts",
                       "06_Cluster_Sizes_all_features/cluster_complete_data", sep="/")

#########################################################################################################
# BRCA
brca_cluster5 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/1_1_1.txt'))))
brca_cluster4 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/1_2_1.txt'))))
brca_cluster3 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/1_3_1.txt'))))
brca_cluster2 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/1_4_1.txt'))))

summary(brca_cluster5)
summary(brca_cluster4)
summary(brca_cluster3)
summary(brca_cluster2)


# KIRC
kirc_cluster5 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/2_1_1.txt'))))
kirc_cluster4 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/2_2_1.txt'))))
kirc_cluster3 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/2_3_1.txt'))))
kirc_cluster2 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/2_4_1.txt'))))

summary(kirc_cluster5)
summary(kirc_cluster4)
summary(kirc_cluster3)
summary(kirc_cluster2)

# KIRC (combat)
combat_cluster5 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/6_5.txt'))))
combat_cluster4 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/6_4.txt'))))
combat_cluster3 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/6_3.txt'))))

summary(combat_cluster5)
summary(combat_cluster4)
summary(combat_cluster3)


# LIHC
lihc_cluster5 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/3_1_1.txt'))))
lihc_cluster4 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/3_2_1.txt'))))
lihc_cluster3 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/3_3_1.txt'))))
lihc_cluster2 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/3_4_1.txt'))))

summary(lihc_cluster5)
summary(lihc_cluster4)
summary(lihc_cluster3)
summary(lihc_cluster2)


# LUAD
luad_cluster5 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/4_1_1.txt'))))
luad_cluster4 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/4_2_1.txt'))))
luad_cluster3 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/4_3_1.txt'))))
luad_cluster2 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/4_4_1.txt'))))

summary(luad_cluster5)
summary(luad_cluster4)
summary(luad_cluster3)
summary(luad_cluster2)


# LUSC
lusc_cluster5 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/5_1_1.txt'))))
lusc_cluster4 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/5_2_1.txt'))))
lusc_cluster3 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/5_3_1.txt'))))
lusc_cluster2 <- as.factor(t(read.table(paste0(testFolder, '/Results/Labels/5_4_1.txt'))))

summary(lusc_cluster5)
summary(lusc_cluster4)
summary(lusc_cluster3)
summary(lusc_cluster2)








