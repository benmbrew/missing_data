#####################################################################################
# this script will look at the final ranked results of intersection and union

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

# initialize fixed variables 
cancer <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")

################################################################################################
# Load data 
int_results <- read.csv(paste0(results_folder, '/missing_data_int.csv'))
union_results <- read.csv(paste0(results_folder, '/missing_data_union.csv'))

################################################################################################
## add rank to for each imputation and clustering method and combination of methods
int_results$X <- NULL

# create a variable for imputation and cluster
createColumn <- 
  
  function(data, impute) {
    
    temp <- strsplit(as.character(data$method), '.', fixed = TRUE)
    if (impute) {
      temp.2 <- lapply(temp, function(x) x[length(x)])
      data$imputation <- temp.2
    } else {
      temp.2 <- lapply(temp, function(x) x[length(x) - 1])
      data$cluster_method <- temp.2
    }
    return(data)
    
  }

int_results <- createColumn(int_results, impute = TRUE)
int_results <- createColumn(int_results, impute = FALSE)

union_results <- createColumn(union_results, impute = TRUE)
union_results <- createColumn(union_results, impute = FALSE)


# get ranking for each method
rankMethods <- function(data, complete) {
  
  if (complete) {
    
    data <- transform(data, 
                      acc_nmi_rank_int = ave(acc_nmi, 
                                             cancer,
                                             cluster,
                                             FUN = function(x) rank(-x, ties.method = "min")),
                      acc_rank_int = ave(acc, 
                                             cancer,
                                             cluster,
                                             FUN = function(x) rank(-x, ties.method = "min")),
                      nmi_rank_int = ave(nmi, 
                                         cancer,
                                         cluster,
                                         FUN = function(x) rank(-x, ties.method = "min")),
                      pval_ci_rank_int = ave(pval_ci, 
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

int_rank <- rankMethods(int_results, complete = TRUE)
union_rank <- rankMethods(union_results, complete = FALSE)

# unlist and make factor for group by later.
int_rank$cluster_method <- as.factor(unlist(int_rank$cluster_method))
int_rank$imputation <- as.factor(unlist(int_rank$imputation))

union_rank$cluster_method <- as.factor(unlist(union_rank$cluster_method))
union_rank$imputation <- as.factor(unlist(union_rank$imputation))

#######################################################################################################
# optimal clusters are 3, 2, 4, 4, 4

#intersection
brca_int <- int_rank[int_rank$cancer == 'BRCA' & int_rank$cluster == 3,]
kirc_int <- int_rank[int_rank$cancer == 'combat' & int_rank$cluster == 2,]
lihc_int <- int_rank[int_rank$cancer == 'LIHC' & int_rank$cluster == 4,]
luad_int <- int_rank[int_rank$cancer == 'LUAD' & int_rank$cluster == 4,]
lusc_int <- int_rank[int_rank$cancer == 'LUSC' & int_rank$cluster == 4,]

int_final <- rbind(brca_int, kirc_int, lihc_int, luad_int, lusc_int)

# union
brca_union <- union_rank[union_rank$cancer == 'BRCA' & union_rank$cluster == 3,]
kirc_union <- union_rank[union_rank$cancer == 'combat' & union_rank$cluster == 2,]
lihc_union <- union_rank[union_rank$cancer == 'LIHC' & union_rank$cluster == 4,]
luad_union <- union_rank[union_rank$cancer == 'LUAD' & union_rank$cluster == 4,]
lusc_union <- union_rank[union_rank$cancer == 'LUSC' & union_rank$cluster == 4,]

union_final <- rbind(brca_union, kirc_union, lihc_union, luad_union, lusc_union)



########################################################################################################
#### METHOD

# group by method, get raning for acc_nmi
top_int_acc <- int_final %>%
  group_by(method) %>%
  summarise(rank1 = sum(acc_rank_int == 1),
            rank2 = sum(acc_rank_int == 2),
            rank3 = sum(acc_rank_int == 3),
            rank4 = sum(acc_rank_int == 4),
            rank5 = sum(acc_rank_int == 5),
            rank6 = sum(acc_rank_int == 6),
            rank7 = sum(acc_rank_int == 7))

# melt top_int 
top_int_acc_melt <- melt(top_int_acc, id.vars = 'method')

# plot 
ggplot(data = top_int_acc_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('Method') + ylab('') + 
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey8", "grey15", "grey30", "grey45", "grey60", "grey75", "grey90"), 
                    name="ACC and NMI Ranking",
                    breaks=c("rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7"),
                    labels=c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7")) + 
  scale_y_continuous(breaks=c(1,2,3,4,5)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="white"), 
        panel.border = element_rect(fill = NA, colour = 'grey50'),
        panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)


#### METHOD

# group by method, get raning for acc_nmi
top_int_acc_nmi <- int_final %>%
  group_by(method) %>%
  summarise(rank1 = sum(acc_nmi_rank_int == 1),
            rank2 = sum(acc_nmi_rank_int == 2),
            rank3 = sum(acc_nmi_rank_int == 3),
            rank4 = sum(acc_nmi_rank_int == 4),
            rank5 = sum(acc_nmi_rank_int == 5),
            rank6 = sum(acc_nmi_rank_int == 6),
            rank7 = sum(acc_nmi_rank_int == 7))

# melt top_int 
top_int_acc_nmi_melt <- melt(top_int_acc_nmi, id.vars = 'method')

# plot 
ggplot(data = top_int_acc_nmi_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('Method') + ylab('') + 
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey8", "grey15", "grey30", "grey45", "grey60", "grey75", "grey90"), 
                    name="ACC and NMI Ranking",
                    breaks=c("rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7"),
                    labels=c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7")) + 
  scale_y_continuous(breaks=c(1,2,3,4,5)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="white"), 
        panel.border = element_rect(fill = NA, colour = 'grey50'),
        panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)

# regress nad median
############################################################################################################
# group by method, get rnaking for total
top_int_total <- int_final %>%
  group_by(method) %>%
  summarise(rank1 = sum(total_rank_int == 1),
            rank2 = sum(total_rank_int == 2),
            rank3 = sum(total_rank_int == 3),
            rank4 = sum(total_rank_int == 4),
            rank5 = sum(total_rank_int == 5),
            rank6 = sum(total_rank_int == 6),
            rank7 = sum(total_rank_int == 7))

# melt top_int 
top_int_total_melt <- melt(top_int_total, id.vars = 'method')

# plot 
ggplot(data = top_int_total_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('Method') + ylab('') + 
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey8", "grey15", "grey30", "grey45", "grey60", "grey75", "grey90"), 
                    name="Total Ranking",
                    breaks=c("rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7"),
                    labels=c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7")) + 
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="white"), 
        panel.border = element_rect(fill = NA, colour = 'grey50'),
        panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)

# regress and median

# group by method, get rnaking for pval ci
top_int_pval_ci <- int_final %>%
  group_by(method) %>%
  summarise(rank1 = sum(pval_ci_rank_int == 1),
            rank2 = sum(pval_ci_rank_int == 2),
            rank3 = sum(pval_ci_rank_int == 3),
            rank4 = sum(pval_ci_rank_int == 4),
            rank5 = sum(pval_ci_rank_int == 5),
            rank6 = sum(pval_ci_rank_int == 6),
            rank7 = sum(pval_ci_rank_int == 7))

# melt top_int 
top_int_pval_ci_melt <- melt(top_int_pval_ci, id.vars = 'method')

# plot 
ggplot(data = top_int_pval_ci_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('Method') + ylab('') + 
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey8", "grey15", "grey30", "grey45", "grey60", "grey75", "grey90"), 
                    name="Pval and CI Ranking",
                    breaks=c("rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7"),
                    labels=c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7")) + 
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="white"), 
        panel.border = element_rect(fill = NA, colour = 'grey50'),
        panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)

# regress, knn, lls

# group by method, get rnaking for pval 
top_int_pval <- int_final %>%
  group_by(method) %>%
  summarise(rank1 = sum(pval_rank_int == 1),
            rank2 = sum(pval_rank_int == 2),
            rank3 = sum(pval_rank_int == 3),
            rank4 = sum(pval_rank_int == 4),
            rank5 = sum(pval_rank_int == 5),
            rank6 = sum(pval_rank_int == 6),
            rank7 = sum(pval_rank_int == 7))

# melt top_int 
top_int_pval_melt <- melt(top_int_pval, id.vars = 'method')

# plot 
ggplot(data = top_int_pval_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('Method') + ylab('') + 
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey8", "grey15", "grey30", "grey45", "grey60", "grey75", "grey90"), 
                    name="Pval Ranking",
                    breaks=c("rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7"),
                    labels=c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7")) + 
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="white"), 
        panel.border = element_rect(fill = NA, colour = 'grey50'),
        panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)


# pval, lls
############################################################################################################

# group by method, get ranking for union
top_union_pval_ci <- union_final %>%
  group_by(method) %>%
  summarise(rank1 = sum(pval_ci_rank_union == 1),
            rank2 = sum(pval_ci_rank_union == 2),
            rank3 = sum(pval_ci_rank_union == 3),
            rank4 = sum(pval_ci_rank_union == 4),
            rank5 = sum(pval_ci_rank_union == 5),
            rank6 = sum(pval_ci_rank_union == 6),
            rank7 = sum(pval_ci_rank_union == 7))

# melt top_int 
top_union_pval_ci_melt <- melt(top_union_pval_ci, id.vars = 'method')

# plot 
ggplot(data = top_union_pval_ci_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('Method') + ylab('') + 
  ggtitle('Union') +
  scale_fill_manual(values = c("grey8", "grey15", "grey30", "grey45", "grey60", "grey75", "grey90"), 
                    name="Pval and CI Ranking on Union",
                    breaks=c("rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7"),
                    labels=c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7")) + 
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="white"), 
        panel.border = element_rect(fill = NA, colour = 'grey50'),
        panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)

# regress and lls

# group by method, get ranking for union pval
top_union_pval <- union_final %>%
  group_by(method) %>%
  summarise(rank1 = sum(pval_rank_union == 1),
            rank2 = sum(pval_rank_union == 2),
            rank3 = sum(pval_rank_union == 3),
            rank4 = sum(pval_rank_union == 4),
            rank5 = sum(pval_rank_union == 5),
            rank6 = sum(pval_rank_union == 6),
            rank7 = sum(pval_rank_union == 7))

# melt top_int 
top_union_pval_melt <- melt(top_union_pval, id.vars = 'method')

# plot 
ggplot(data = top_union_pval_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.8)  +  xlab('Method') + ylab('') + 
  ggtitle('Union') +
  scale_fill_manual(values = c("grey8", "grey15", "grey30", "grey45", "grey60", "grey75", "grey90"), 
                    name="Pval Ranking on Union",
                    breaks=c("rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7"),
                    labels=c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Rank 7")) + 
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="white"), 
        panel.border = element_rect(fill = NA, colour = 'grey50'),
        panel.grid.major=element_line(colour="grey50",size=0.2, linetype = 'dashed'), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) + geom_hline(yintercept=0)

##################################################################################################
# Load in raw

