# Return the imputed data and performance metrics for the method.
# impute: function used to impute the incomplete data
# incompleteData: list of matrices with missing values
# removedData: list of matrices which contain the values missing in
# incompleteData

evaluateImputation <- function(impute, incompleteData, removedData) {
  # Impute the incomplete data
  runTime <- system.time(imputedData <- impute(incompleteData))[3]
  
  # Evaluate the imputed data
  # Note: the metrics for each data type is calculated separately
  nrmse <- sapply(1:length(incompleteData), function(i)
    imputationErrors(removedData[[i]], imputedData[[i]])$nrmse)
  
  results <- c(nrmse, runTime)

  return(list(imputedData=imputedData, results=results))
}