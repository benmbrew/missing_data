# Return metrics to evaluate clustering the incomplete data using
# similarity methods.
# similarity: function used to impute missing affinity values
# data: list of matrices containing missing values
# completeLabels: numeric vector of labels derived from clustering the
# complete data
# clinicalData: data frame of TCGA clinical data
# sampleRows: TRUE iff samples are contained in the rows of the data

evaluateSimilarityCustom <- function(similarity, data, completeLabels,
                               clinicalData, intersectInd,
                               sampleRows) {
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
  
  # Replace the missing similarities
  runTime <- system.time(affinities <- similarity(affinities))[3]
  
  # Fuse the affinity matrices
  fusedMatrix <- SNF(affinities)
  
  # Cluster the fused matrix
  # get number of clusters: This function estimates the number of clusters given 
  # the two huristics given in the supplementary materials of our nature method paper W 
  # is the similarity graph NUMC is a vector which contains the possible choices of number of clusters.
  numClusEstimates <- unlist(estimateNumberOfClustersGivenGraph(fusedMatrix, 2:10))
  numClus <- max(numClusEstimates[c(1,3)])
  # chooses 2
  labels <- spectralClustering(fusedMatrix, numClus)
  
  # Evaluate the new labels
  acc <- calAcc(completeLabels, labels)
  nmi <- calNMI(completeLabels, labels)
  surv <- evaluateSurvival(clinicalData, labels)
  
  # Evaluate the new labels on individuals in the intersection
  intLabels <- labels[intersectInd]
  intCompleteLabels <- completeLabels[intersectInd]
  intClinicalData <- clinicalData[intersectInd, ]
  intAcc <- calAcc(intCompleteLabels, intLabels)
  intNmi <- calNMI(intCompleteLabels, intLabels)
  if (length(unique(intLabels))==1) {
    intSurv <- c(NA, NA)
  } else {
    intSurv <- evaluateSurvival(intClinicalData, intLabels)
  }
  
  results <- c(acc, nmi, surv, intAcc, intNmi, intSurv, runTime, numClus)
  
  return(results)
}