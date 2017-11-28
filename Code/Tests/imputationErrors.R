# Return a list of imputation errors.
# removed: a matrix containing the values of entries which will be
# used to calculate the error; all other entries are NA
# imputed: a matrix containing imputed values at the same positions
# as values present in the removed matrix 

imputationErrors <- function(removed, imputed) {
  mse <- 0
  rmse <- 0
  nrmse <- 0
  rae <- 0
  if (!all(is.na(removed))) {
    mse <- mean((removed-imputed)^2, na.rm=TRUE) # MSE
    rmse <- mse^0.5 # RMSE
    nrmse <- rmse/sd(as.vector(removed), na.rm=TRUE) # NRMSE
    
    absError <- abs(removed - imputed)
    relAbsError <- absError / abs(removed)
    rae <- mean(relAbsError, na.rm=TRUE) # RAE
  }
  
  return(list(mse=mse, rmse=rmse, nrmse=nrmse, rae=rae))
}