
##########
# load libraries
##########
library(tidyverse)

method_length = 15

##########
# load data 
##########
intersection <- read.csv('~/Desktop/final_results_intersection.csv')
union <- read.csv('~/Desktop/final_results_union.csv')

##########
# keep only necessary columns
##########
intersection <- intersection[, c('cancer', 'method', 'acc', 'nmi', 'pval', 'ci')]

intersection[, 3:ncol(intersection)] <- method_length - intersection[, 3:ncol(intersection)]

ranked_int <- intersection %>% 
  group_by(cancer) %>%
  mutate(rank_acc = floor(rank(acc)),
         rank_nmi = floor(rank(nmi)),
         rank_pval = floor(rank(pval)),
         rank_ci = floor(rank(ci)))


write.csv(ranked_int, '~/Desktop/int_rank.csv')

##########
# keep only necessary columns
##########
union <- union[, c('cancer', 'method', 'pval', 'ci')]

union[, 3:ncol(union)] <- method_length - union[, 3:ncol(union)]

ranked_union <- union %>% 
  group_by(cancer) %>%
  mutate(rank_pval = floor(rank(pval)),
         rank_ci = floor(rank(ci)))

write.csv(ranked_union, '~/Desktop/union_rank.csv')




