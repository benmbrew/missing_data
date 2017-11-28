lsaImputation <- function(incompleteData, sampleRows) {
  # Variables which have not been initialized:
  # imputedFile, incompleteFile, projectFolder, jvmGBLimit
  
  # Transpose the data for LSA
  # LSA requires samples in rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    incompleteData <- lapply(incompleteData, t)
  }
  
  # Concatenate the data by columns because the samples are in rows
  # Row and column names are required for the LSA code
  concatenatedIncompleteData <- do.call(cbind, incompleteData) # cols
  concatenatedRownames <- rownames(concatenatedIncompleteData) # rows
  colnames(concatenatedIncompleteData) <-
    as.character(1:ncol(concatenatedIncompleteData))
  rownames(concatenatedIncompleteData) <-
    as.character(1:nrow(concatenatedIncompleteData))
  
  # Save the incomplete data to a file
  cat("xxx","\t", file=incompleteFile, sep="")
  write.table(concatenatedIncompleteData, incompleteFile,
              col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE,
              na="NULL",append=TRUE)
  
  # Impute the missing data
  lsaCodePath <- paste(projectFolder, "Code/Impute/LSimpute.jar",
                       sep="/")
  jvmMemLimit <- paste("-Xmx", jvmGBLimit, "G", sep="")
  javaCommand <- paste("java -jar", jvmMemLimit, "-server ",
                       lsaCodePath, incompleteFile, imputedFile, "6",
                       sep=" ")
  system(javaCommand)
  
  # Read the imputed file
  concatenatedImputedData <-
    data.matrix(read.table(imputedFile, header=TRUE, sep="\t"))
  concatenatedImputedData <- concatenatedImputedData[, -1]
  rownames(concatenatedImputedData) <- concatenatedRownames # rows
  
  # Remove the files used for imputation
  system(paste("rm", imputedFile, sep=" "))
  system(paste("rm", incompleteFile, sep=" "))
  
  # Split the concatenated data by columns
  featureSize <- sapply(incompleteData, ncol) # cols
  imputedData <- splitConcatenatedData(concatenatedImputedData, 
                                       featureSize,
                                       sampleRowsRequired)
  
  # Reorient the data to its original form
  if (transposeData) {
    imputedData <- lapply(imputedData, t)
  }
  
  return(imputedData)
}