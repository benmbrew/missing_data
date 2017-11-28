# Return the labels from SNF clustering, given the list of data and
# the number of clusters.

SNFClusteringType <- function(data, numClus, sampleRows) {
  # Transpose data for the dist function
  # The dist function takes distances between rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- t(data)
  }
  
  # Calculate the distance between samples
  distances <- as.matrix(dist(data))
  
  # Convert the distances to affinities
  affinities <- affinityMatrix(distances)
  
  # Fuse the affinity matrices
  #fusedMatrix <- SNF(affinities)
  
  # Cluster the fused matrix
  labels <- spectralClustering(affinities, numClus)
  
  return(labels)
}