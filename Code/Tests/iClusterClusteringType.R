# Return the labels from iCluster clustering, given the list of data
# and the number of clusters.

iClusterClusteringType <- function(data, numClus, sampleRows) {
  # Transpose the data for the iCluster function.
  # The iCluster function requires samples in rows.
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- t(data)
  }
  
  # Fit the model using a range of lambda parameters
  # K = number of eigen-features = numClus - 1
  type <- rep.int("gaussian", 3)
  n.lambda <- 35
  scale.lambda <- rep.int(1, 3)
  K <- numClus - 1
  cvFit <- tune.iClusterPlus(cpus=numCores, data, type=type, K=K,
                             n.lambda=n.lambda,
                             scale.lambda=scale.lambda)
  
  # Cluster the data using lambda which minimizes the BIC
  BIC <- getBIC(list(cvFit))
  minBICind <- which.min(BIC)
  labels <- cvFit$fit[[minBICind]]$clusters
  
  # Return the resulting labels
  return(labels)
}