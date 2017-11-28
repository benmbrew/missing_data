# Return a list of matrices derived from splitting the features of a
# concatenated matrix.
# concatenatedData: a matrix which will be split
# featureSize: a vector containing the number of features contained in
# the resulting matrices
# sampleRows: TRUE iff samples are in the rows of the data

splitConcatenatedData <- function(concatenatedData, featureSize,
                                  sampleRows) {
  # Transpose the data so that the samples are contained in the rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    concatenatedData <- t(concatenatedData)
  }
  
  # Split the data by columns because the samples are in the rows
  splitSize <- length(featureSize)
  splitData <- vector("list", splitSize)
  for (v in 1:splitSize) {
    startInd <- sum(featureSize[0:(v-1)]) + 1
    endInd <- sum(featureSize[1:v])
    splitData[[v]] <- concatenatedData[, startInd:endInd]
  }
  
  # Reorient the data to its original form
  if (transposeData) {
    splitData <- lapply(splitData, t)
  }
  
  return(splitData)
}