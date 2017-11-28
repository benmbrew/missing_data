randomImputation <- function(incompleteData, sampleRows, seed=0) {
  # Transpose the data for the method
  # Method requires samples in the columns
  sampleRowsRequired <- FALSE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    incompleteData <- lapply(incompleteData, t)
  }
  
  # Impute missing values in a vector using randomly selected
  # non-missing values
  randomVectorImputation <- function(incompleteVector) {
    imputedVector <- incompleteVector
    missingInd <- is.na(incompleteVector)
    imputeInd <- sample(which(!missingInd), sum(missingInd),
                        replace=TRUE)
    imputedVector[missingInd] <- incompleteVector[imputeInd]
    return(imputedVector)
  }
  
  # Impute missing values in a matrix using randomly selected
  # non-missing values from the same row
  randomMatrixImputation <- function(incompleteMatrix) {
    t(apply(incompleteMatrix, 1, randomVectorImputation))
  }
  
  # Impute the missing data
  imputedData <- lapply(incompleteData, randomMatrixImputation)
  
  # Reorient the data to its original form
  if (transposeData) {
    imputedData <- lapply(imputedData, t)
  }
  
  return(imputedData)
}