# Return the labels from hierarchical clustering the concatenated
# data, given the list of data and the number of clusters.
# Distance metric: 1 - correlation
# Linkage function: average distance

hierarchicalClusteringType <- function(data, numClus, sampleRows) {
  # Transpose the data for the cor function
  # The cor function takes correlations between columns
  sampleRowsRequired <- FALSE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- t(data)
  }
  
  # Concatenate the data which is aligned by column
  concatenatedData <- data
  
  # Calculate the correlation between samples
  # The cor function takes correlations between columns
  correlation <- cor(concatenatedData)
  
  # Convert the correlation to a distance object
  distance <- as.dist(1 - correlation)
  
  # Perform clustering
  hclustFit <- hclust(distance, method="average")
  labels <- cutree(hclustFit, k=numClus)
  
  # Return the resulting labels
  return(labels)
}