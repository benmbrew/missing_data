# Return the labels from SNF clustering, given the list of data and
# the number of clusters.

SNFClustering <- function(data, numClus, sampleRows) {
  # Transpose data for the dist function
  # The dist function takes distances between rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- lapply(data, t)
  }
  
  # Calculate the distance between samples
  distances <- lapply(data, function(x) as.matrix(dist(x)))
  
  # Convert the distances to affinities
  affinities <- lapply(distances, affinityMatrix)
  temp <- affinities[[1]]
  # Fuse the affinity matrices
  fusedMatrix <- SNF(affinities)
  
  # Cluster the fused matrix
  labels <- spectralClustering(fusedMatrix, numClus)
  
  return(labels)
}