regressionSimilarity <- function(Wall) {
  nViews <- length(Wall)
  newWall <- vector("list", nViews)
  
  for (v in 1:nViews) {
    W <- Wall[[v]]
    
    if (any(is.na(W))) {
      # Replace the diagonal missing values.
      y <- diag(W)
      x <- Reduce(cbind, lapply(Wall[-v], diag))
      x <- as.matrix(x)
      # Replace the missing values in the features with the mean similarity.
      for (k in 1:(nViews-1)) {
        x[is.na(x[, k]), k] <- mean(x[, k], na.rm=T)
      }
      diagonalData <- as.data.frame(cbind(y, x))
      diagonalModel <- lm(y~., diagonalData)
      diag(W)[is.na(y)] <- predict(diagonalModel, diagonalData[is.na(y), ])
      
      # Replace the off-diagonal missing values.
      offDiag <- function(W) {
        diag(W) <- Inf
        W <- W[!is.infinite(W)]
        return(W)
      }
      y <- offDiag(W)
      x <- Reduce(cbind, lapply(Wall[-v], offDiag))
      x <- as.matrix(x)
      for (k in 1:(nViews-1)) {
        x[is.na(x[, k]), k] <- mean(x[, k], na.rm=T)
      }
      offDiagonalData <- as.data.frame(cbind(y, x))
      offDiagonalModel <- lm(y~., offDiagonalData)
      missingInd <- is.na(W)
      W[missingInd] <- predict(offDiagonalModel, offDiagonalData[is.na(y), ])
    }
    
    newWall[[v]] <- W
  }
  
  return(newWall)
}