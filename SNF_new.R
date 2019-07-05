load('~/Desktop/new_snf.RData')


similarity
data
completeLabels <- completeLabels[dataInd]
clinicalData <- clinicalData[dataInd, ]
sampleRows <- F
evaluateSimilarity <- function(similarity, data, completeLabels,
                               clinicalData, sampleRows) {
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
  # SNF algorithm creating NAs
  fusedMatrix <- SNFtool::SNF(affinities)
  
  # Cluster the fused matrix
  numClus <- length(unique(completeLabels))
  labels <- spectralClustering(fusedMatrix, numClus)
  
  # Evaluate the new labels
  acc <- calAcc(completeLabels, labels)
  nmi <- calNMI(completeLabels, labels)
  surv <- evaluateSurvival(clinicalData, labels)
  
  #   # Evaluate the new labels on individuals in the intersection
  #   intLabels <- labels[intersectInd]
  #   intCompleteLabels <- completeLabels[intersectInd]
  #   intClinicalData <- clinicalData[intersectInd, ]
  #   intAcc <- calAcc(intCompleteLabels, intLabels)
  #   intNmi <- calNMI(intCompleteLabels, intLabels)
  #   if (length(unique(intLabels))==1) {
  #     intSurv <- c(NA, NA)
  #   } else {
  #     intSurv <- evaluateSurvival(intClinicalData, intLabels)
  #   }
  
  results <- c(acc, nmi, surv) #intAcc, intNmi, intSurv, runTime)
  
  return(results)
}

saveRDS(affinities, '~/Desktop/affinities.rda')

Wall <- affinities

length(affinities)

function (Wall, K = 20, t = 20) 
{
  LW = length(Wall)
  normalize <- function(X) {
    X <- X/(2 * (rowSums(X) - diag(X)))
    diag(X) <- 0.5
    return(X)
  }
  # HERE is problem
  newW <- vector("list", LW)
  nextW <- vector("list", LW)
  for (i in 1:LW) {
    # code
    Wall[[i]] = normalize(Wall[[i]])
    Wall[[i]] = (Wall[[i]] + t(Wall[[i]]))/2
  }
  for (i in 1:LW) {
    newW[[i]] = (.dominateset(Wall[[i]], K))
  }
  for (i in 1:t) {
    for (j in 1:LW) {
      sumWJ = matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
      for (k in 1:LW) {
        if (k != j) {
          sumWJ = sumWJ + Wall[[k]]
        }
      }
      nextW[[j]] = newW[[j]] %*% (sumWJ/(LW - 1)) %*% t(newW[[j]])
    }
    for (j in 1:LW) {
      Wall[[j]] <- normalize(nextW[[j]])
      Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2
    }
  }
  W = matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
  for (i in 1:LW) {
    W = W + Wall[[i]]
  }
  W = W/LW
  W = normalize(W)
  W = (W + t(W))/2
  return(W)
}