llsImputation <- function(incompleteData, sampleRows) {
  # Transpose the data for LLS
  # LLS requires samples in the rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    incompleteData <- lapply(incompleteData, t)
  }
  
  # Concatenate the data by columns because the samples are in rows
  concatenatedIncompleteData <- do.call(cbind, incompleteData)
  
  # Impute the missing data
  concatenatedImputedData <- llsImpute(concatenatedIncompleteData,
                                       allVariables=TRUE)
  
  # Split the concatenated data by columns
  featureSize <- sapply(incompleteData, ncol)
  imputedData <- splitConcatenatedData(concatenatedImputedData,
                                       featureSize,
                                       sampleRowsRequired)
  
  # Reorient the data to its original form
  if (transposeData) {
    imputedData <- lapply(imputedData, t)
  }
  
  return(imputedData)
}