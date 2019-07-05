# Load libraries
library(SNFtool)

#  create a list to store data
data <- list()

# Load the biological data
data[[1]] <- as.matrix(read.delim('Data/habib/GLOM_98_GE_ModAv.txt'))
data[[2]] <- as.matrix(read.delim('Data/habib/TUB_98_GE_ModAv.txt'))

# subset data to run locally
data[[1]] <- data[[1]][1:100, 1:30]
data[[2]] <- data[[2]][1:100, 1:30]

sampleRows <- FALSE
# Transpose data for the dist function
# The dist function takes distances between rows
sampleRowsRequired <- TRUE
transposeData <- sampleRows != sampleRowsRequired
if (transposeData) {
  data <- lapply(data, t)
}

# store the first row (features) and remove
features <- list()
features[[1]] <- data[[1]][1,]
features[[2]] <- data[[2]][1,]
data[[1]] <- data[[1]][-1,]
data[[2]] <- data[[2]][-1,]

# convert to numeric
temp_1 <- data[[1]]
temp_1 <- apply(temp_1, 2, function(x) as.numeric(x))
temp_2 <- data[[2]]
temp_2 <- apply(temp_2, 2, function(x) as.numeric(x))

# store back in lisst
data[[1]] <- temp_1
data[[2]] <- temp_2

rm(temp_1, temp_2)
# Calculate the distance between samples
distances <- lapply(data, function(x) as.matrix(dist(x)))

# Convert the distances to affinities
affinities <- lapply(distances, affinityMatrix)
# Fuse the affinity matrices
fusedMatrix <- SNF(affinities)

#save.image('~/Desktop/temp_snf.RData')
# load('~/Desktop/temp_snf.RData')
# apply nmf clustering

# there are a lot of different ways to specify the desired rank.
nmf_dat <- NMF::nmf(fusedMatrix, rank = 3)

# W matrix
w <- nmf_dat@fit@W

# V matrix
h<- nmf_dat$fit@H

# not sure how to get clusters from this?
v < w%*%h




