# Given a list of matrices, L, with column names, return a list of matrices
# which are the same as those in L, except they contain additional columns of
# NAs for each column name which only appears in other matrices in L.  The
# columns of the new matrices are sorted in alphabetical order, by their column
# name.

columnUnion <- function(L) {
  n <- length(L)
  cnames <- colnames(L[[1]])
  
  for (i in 1:n) {
    cnames <- union(cnames, colnames(L[[i]]))
  }
  cnames <- sort(cnames)

  unionL <- vector("list", n)
  for (i in 1:n) {
    unionL[[i]] <- L[[i]][, match(cnames, colnames(L[[i]])), drop=FALSE]
    colnames(unionL[[i]]) <- cnames
  }
  unionL
}

sampleUnion <- function(data, sampleRows) {
  if (sampleRows) {
    data <- lapply(data, t)
  }
  
  data <- columnUnion(data)
  
  if (sampleRows) {
    data <- lapply(data, t)
  }
  
  return(data)
}