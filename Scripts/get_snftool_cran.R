## ======================================================================
## Step 1: Download all log files
## ======================================================================
library(data.table)

# Here's an easy way to get all the URLs in R
start <- as.Date('2016-01-01')
today <- as.Date('2016-06-01')

all_days <- seq(start, today, by = 'day')

year <- as.POSIXlt(all_days)$year + 1900
urls <- paste0('http://cran-logs.rstudio.com/', year, '/', all_days, '.csv.gz')

# only download the files you don't have:
missing_days <- setdiff(as.character(all_days), tools::file_path_sans_ext(dir("CRANlogs"), TRUE))

#dir.create("CRANlogs")
for (i in 1:length(missing_days)) {
  print(paste0(i, "/", length(missing_days)))
  download.file(urls[i], paste0('CRANlogs/', missing_days[i], '.csv.gz'))
}


## ======================================================================
## Step 2: Load single data files into one big data.table
## ======================================================================

file_list <- list.files("CRANlogs", full.names=TRUE)


logs <- list()
for (file in file_list) {
  print(paste("Reading", file, "..."))
  logs[[file]] <- read.table(file, header = TRUE, sep = ",", quote = "\"",
                             dec = ".", fill = TRUE, comment.char = "", as.is=TRUE)
}

# rbind together all files
dat <- rbindlist(logs)

# add some keys and define variable types
dat[, date:=as.Date(date)]
dat[, package:=factor(package)]
dat[, country:=factor(country)]
dat[, weekday:=weekdays(date)]
dat[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]

setkey(dat, package, date, week, country)
rm(logs)

save(dat, file="CRANlogs/second.RData")

# for later analyses: load the saved data.table
# load("/home/benbrew/Desktop/rdata_cran/CRANlogs.RData")
# load("/home/benbrew/Desktop/rdata_cran/CRANlogs1.RData")
# load("/home/benbrew/Desktop/rdata_cran/CRANlogs2.RData")
# load("/home/benbrew/Desktop/rdata_cran/CRANlogs3.RData")



# save.image('/home/benbrew/Desktop/data_frame.RData')
# =====================================================================
  ## Step 3: Analyze it!
  ## ======================================================================


library(ggplot2)
library(plyr)

str(dat)

# Overall downloads of packages
d4 <- dat[, length(week), by=package]
d4 <- d4[order(V1), ]
d4[package=="SNFtool", ]

# # plot 1: Compare downloads of selected packages on a weekly basis
# agg1 <- data[J(c("SNFtool")), length(unique(ip_id)), by=c("week", "package")]
# 
# ggplot(agg1, aes(x=week, y=V1, color=package, group=package)) + 
#   geom_line() + ylab("Downloads") + theme_bw() + 
#   theme(axis.text.x  = element_text(angle=90, size=8, vjust=0.5))
# 

load('/home/benbrew/Desktop/data_frame.RData')

# convert to data frame
d1 <- as.data.frame(d1)
d2 <- as.data.frame(d2)
d3 <- as.data.frame(d3)
d4 <- as.data.frame(d4)

# add time to d1-d4
d1$time <- '2014'
d2$time <- '2015'
d3$time <- '2015'
d4$time <- '2016'


# rbind 

data <- rbind(d1, d2, d3, d4)
data <- data[data$package == 'SNFtool',]
data[2, 2] <- data[2,2] + data[3,2]
data <- data[-3,]


# remove NAs
data <- data[complete.cases(data),]

ggplot(data = data, aes(time, V1)) + geom_bar(stat = 'identity', fill = 'grey40', alpha = 0.6) + 
  xlab('') + ylab('# of Downloads') + ggtitle('SNFtool Downloads') + 
  geom_text(aes(label=V1), vjust=1.6, color="white", size=7) + 
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major = element_line(linetype = 'dashed', colour = "#F0F0F0"),
        plot.background=element_rect(fill="#F0F0F0"),
        panel.grid.major=element_line(colour="#F0F0F0",size=.75), axis.ticks=element_blank(),
        legend.position="right", legend.title = element_blank(),
        legend.background = element_rect(fill="#F0F0F0"),
        plot.title=element_text(face="bold",hjust=0,vjust=2,colour="#535353",size=25),
        axis.text.x=element_text(size=10,colour="#535353",face="bold", angle = 0, hjust = 0.5),
        axis.text.y=element_text(size=15,colour="#535353",face="bold"),
        axis.title.y=element_text(size=15,colour="#535353",face="bold",vjust=1.5),
        axis.title.x=element_text(size=15,colour="#535353",face="bold",vjust=-.5),
        plot.margin = unit(c(1, 1, .5, .7), "cm"))
