## ======================================================================
## Step 1: Download all log files
## ======================================================================
# setwd('/home/benbrew')
# # Here's an easy way to get all the URLs in R
# start <- as.Date('2016-01-01')
# today <- as.Date('2016-06-20')
# 
# all_days <- seq(start, today, by = 'day')
# 
# year <- as.POSIXlt(all_days)$year + 1900
# urls <- paste0('http://cran-logs.rstudio.com/', year, '/', all_days, '.csv.gz')
# 
# # only download the files you don't have:
# missing_days <- setdiff(as.character(all_days), tools::file_path_sans_ext(dir("CRANlogs"), TRUE))
# 
# #dir.create("CRANlogs")
# for (i in 1:length(missing_days)) {
#   print(paste0(i, "/", length(missing_days)))
#   download.file(urls[i], paste0('CRANlogs/', missing_days[i], '.csv.gz'))
# }
# 
# 
# ## ======================================================================
# ## Step 2: Load single data files into one big data.table
# ## ======================================================================
# 
# file_list <- list.files("CRANlogs", full.names=TRUE)
# 
# 
# logs <- list()
# for (file in file_list) {
#   print(paste("Reading", file, "..."))
#   logs[[file]] <- read.table(file, header = TRUE, sep = ",", quote = "\"",
#                              dec = ".", fill = TRUE, comment.char = "", as.is=TRUE)
# }
# 
# # rbind together all files
# dat <- rbindlist(logs)
# 
# # add some keys and define variable types
# dat[, date:=as.Date(date)]
# dat[, package:=factor(package)]
# dat[, country:=factor(country)]
# dat[, weekday:=weekdays(date)]
# dat[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]
# 
# setkey(dat, package, date, week, country)
# rm(logs)
# 
# save(dat, file="CRANlogs/third.RData")

#############################################################################################
# 
# load("/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/SNF/NM_2015/Scripts/06_Results/first.RData")
# # Overall downloads of packages
# dat <- dat[, length(week), by=list(date, package)]
# dat <- dat[package=="SNFtool", ]
# dat1 <- dat
# rm(dat)
# 
# load("/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/SNF/NM_2015/Scripts/06_Results/second.RData")
# # Overall downloads of packages
# dat <- dat[, length(week), by=list(date, package)]
# dat <- dat[package=="SNFtool", ]
# dat2 <- dat
# rm(dat)
# 
# load("/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/SNF/NM_2015/Scripts/06_Results/third.RData")
# # Overall downloads of packages
# dat <- dat[, length(week), by=list(date, package)]
# dat <- dat[package=="SNFtool", ]
# dat3 <- dat
# rm(dat)
# 
# ########################################################################################
# # rbind all three and analyze 
# dat <- as.data.frame(rbind(dat1,dat2,dat3))
# save(dat, 
#      file="/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/SNF/NM_2015/Scripts/06_Results/snf_tool.RData")

load("/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/SNF/NM_2015/Scripts/06_Results/snf_tool.RData")

# Load in images and create new names 
library(ggplot2)
library(plyr)
library(data.table)

# line plot
ggplot(dat, aes(date, V1)) + geom_line(size = 0.6) +
  xlab('Date') + ylab('# of Downloads') + ylim(c(0, 100))

# make cumulative 
dat <- dat[order(dat$date),]
dat$cumulative <- cumsum(dat$V1)

# line plot
ggplot(dat, aes(date, cumulative)) + geom_line(size = 0.6) +
  xlab('Date') + ylab('# of Downloads') 

# create months and years variables 



source("https://bioconductor.org/biocLite.R")
biocLite("minfi")
