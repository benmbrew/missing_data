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
int_results <- read.csv(paste0(results_folder, '/final_results_intersection.csv'))
union_results <- read.csv(paste0(results_folder, '/final_results_union.csv'))

################################################################################################
# Get rid of unnecessary columns 
int_results$cluster <- NULL
int_results$union_rank <- NULL
union_results$clusters <- NULL

################################################################################################
## add rank to for each imputation and clustering method and combination of methods

# create a variable for imputation and cluster
createColumn <- 
  function(data, impute) {
  
  temp <- strsplit(as.character(data$method), '.', fixed = TRUE)
  if (impute) {
    temp.2 <- lapply(temp, function(x) x[length(x)])
    data$imputation <- temp.2
  } else {
    temp.2 <- lapply(temp, function(x) x[length(x) - 1])
    data$cluster <- temp.2
  }
  return(data)
  
}

int_results <- createColumn(int_results, impute = TRUE)
int_results <- createColumn(int_results, impute = FALSE)

union_results <- createColumn(union_results, impute = TRUE)
union_results <- createColumn(union_results, impute = FALSE)


# create ranking for imputation, clustering, and the combination.
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

int_results <- rankMethods(int_results, complete = TRUE)
union_results <- rankMethods(union_results, complete = FALSE)

# unlist and make factor for group by later.
int_results$cluster <- as.factor(unlist(int_results$cluster))
int_results$imputation <- as.factor(unlist(int_results$imputation))

union_results$cluster <- as.factor(unlist(union_results$cluster))
union_results$imputation <- as.factor(unlist(union_results$imputation))


########################################################################################################
#### METHOD
## get top 1, 3, and 5 indicators for intersection and union 

# group by method, get top 1, top 3, top 5 for acc_nmi
top_int_acc <- int_results %>%
  group_by(method) %>%
  summarise(top1_acc_nmi = sum(acc_nmi_rank == 1),
            top3_acc_nmi = sum(acc_nmi_rank < 4),
            top5_acc_nmi = sum(acc_nmi_rank < 6))

# melt top_int 
top_int_acc_melt <- melt(top_int_acc, id.vars = 'method')

# plot 
ggplot(data = top_int_acc_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.7)  +  xlab('Method') + ylab('') + 
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey10", "grey50", "grey90"), name="ACC and NMI Ranking",
                              breaks=c("top1_acc_nmi", "top3_acc_nmi", "top5_acc_nmi"),
                               labels=c("Rank 1", "Top 3", "Top 5")) + scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

############################################################################################################
# group by method, get top 1, top 3, top 5 for pval
top_int_pval <- int_results %>%
  group_by(method) %>%
  summarise(top1_pval = sum(pval_ci_rank == 1),
            top3_pval = sum(pval_ci_rank < 4),
            top5_pval = sum(pval_ci_rank < 6))

# melt top_int 
top_int_pval_melt <- melt(top_int_pval, id.vars = 'method')

# plot 
ggplot(data = top_int_pval_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity')  +  xlab('Method') + ylab('') + 
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey10", "grey50", "grey90"), name="Pval and CI Ranking ",
                    breaks=c("top1_acc_nmi", "top3_acc_nmi", "top5_acc_nmi"),
                    labels=c("Rank 1", "Top 3", "Top 5")) + scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

#############################################################################################################
# group by method, get top 1, top 3, top 5 for pval
top_union <- union_results %>%
  group_by(method) %>%
  summarise(top1_pval = sum(total_rank == 1),
            top3_pval = sum(total_rank < 4),
            top5_pval = sum(total_rank < 6))

# melt top_union 
top_union_melt <- melt(top_union, id.vars = 'method')

# plot 
ggplot(data = top_union_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity')  +  xlab('Method') + ylab('') + scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  ggtitle('Union') +
  scale_fill_manual(values = c("grey10", "grey50", "grey90"), name="Pval and CI Ranking",
                    breaks=c("top1_pval", "top3_pval", "top5_pval"),
                    labels=c("Rank 1", "Top 3", "Top 5")) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

#########################################################################################################
# group by method, get top 1, top 3, top 5 acc_nmi
top_int_acc <- int_results %>%
  group_by(method) %>%
  summarise(top20 = sum(acc_nmi_rank > 0 & acc_nmi_rank <= 3),
            top20_40 = sum(acc_nmi_rank > 3 & acc_nmi_rank <= 6),
            top40_60 = sum(acc_nmi_rank > 6 & acc_nmi_rank <= 9),
            top60_80 = sum(acc_nmi_rank > 9 & acc_nmi_rank <= 12),
            top80_100 = sum(acc_nmi_rank > 12 & acc_nmi_rank <= 15))

# melt top_int 
top_int_acc_melt <- melt(top_int_acc, id.vars = 'method')

# plot 
ggplot(data = top_int_acc_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.7)  +  xlab('Method') + ylab('') + scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey10", "grey30", "grey50", "grey70", "grey90"), name="ACC and NMI Ranking",
                    breaks=c("top20", "top20_40", "top40_60", "top60_80", "top80_100" ),
                    labels=c("1st Quintile", "2nd Quintile", "3rd Quintile", "4th Quintile", "5th Quintile")) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

###########################################################################################################
# group by method, get top 1, top 3, top 5 pval
top_int_pval <- int_results %>%
  group_by(method) %>%
  summarise(top20 = sum(pval_ci_rank > 0 & pval_ci_rank <= 3),
            top20_40 = sum(pval_ci_rank > 3 & pval_ci_rank <= 6),
            top40_60 = sum(pval_ci_rank > 6 & pval_ci_rank <= 9),
            top60_80 = sum(pval_ci_rank > 9 & pval_ci_rank <= 12),
            top80_100 = sum(pval_ci_rank > 12 & pval_ci_rank <= 15))

# melt top_int 
top_int_pval_melt <- melt(top_int_pval, id.vars = 'method')

# plot 
ggplot(data = top_int_pval_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.7)  +  xlab('Method') + ylab('') + scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  ggtitle('Intersection') +
  scale_fill_manual(values = c("grey10", "grey30", "grey50", "grey70", "grey90"), name="Pval and CI Ranking",
                    breaks=c("top20", "top20_40", "top40_60", "top60_80", "top80_100" ),
                    labels=c("1st Quintile", "2nd Quintile", "3rd Quintile", "4th Quintile", "5th Quintile")) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

#######################################################################################################
# group by method, get top 1, top 3, top 5 for pval 
top_union <- union_results %>%
  group_by(method) %>%
  summarise(top20 = sum(total_rank > 0 & total_rank <= 3),
            top20_40 = sum(total_rank > 3 & total_rank <= 6),
            top40_60 = sum(total_rank > 6 & total_rank <= 9),
            top60_80 = sum(total_rank > 9 & total_rank <= 12),
            top80_100 = sum(total_rank > 12 & total_rank <= 15))

# melt top_int 
top_union_melt <- melt(top_union, id.vars = 'method')

# plot 
ggplot(data = top_union_melt, aes(reorder(method, -value), value, fill = variable)) +
  geom_bar(stat = 'identity', alpha = 0.7)  +  xlab('Method') + ylab('') + 
  ggtitle('Union') +
  scale_fill_manual(values = c("grey10", "grey30", "grey50", "grey70", "grey90"), name="Pval and CI Ranking",
                    breaks=c("top20", "top20_40", "top40_60", "top60_80", "top80_100" ),
                    labels=c("1st Quintile", "2nd Quintile", "3rd Quintile", "4th Quintile", "5th Quintile")) +
  theme(panel.background=element_rect(fill="white"), 
        plot.background=element_rect(fill="#F0F0F0"), 
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", #legend.title = element_blank(), 
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=15),
        axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
        axis.text.y=element_text(size=11,colour="#535353",face="bold"),
        axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm")) 

