# Keeping track of the splitting conditions would be easier if there
# was a single variable which could be used to keep track of the logic

knnImputation <- function(incompleteData, sampleRows) {
  # Transpose the data for KNN
  # KNN requires samples in the columns
  sampleRowsRequired <- FALSE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    incompleteData <- lapply(incompleteData, t)
  }
  
  # Concatenate the data by rows because the samples are in columns
  concatenatedIncompleteData <- do.call(rbind, incompleteData)
  
  # Impute the missing data
  concatenatedImputedData <- impute.knn(concatenatedIncompleteData, 
                                        rowmax=1, colmax=1)$data
  
  # Split the concatenated data by rows
  featureSize <- sapply(incompleteData, nrow)
  imputedData <- splitConcatenatedData(concatenatedImputedData,
                                       featureSize,
                                       sampleRowsRequired)
  
  # Reorient the data to its original form
  if (transposeData) {
    imputedData <- lapply(imputedData, t)
  }
  
  return(imputedData)
}