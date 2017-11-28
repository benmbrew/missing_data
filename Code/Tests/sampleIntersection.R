# Given a list of matrices, L, with column names, return a list of matrices
# which are the same as those in L, except they only contain columns whose
# names appear in all of the matrices in L.  The columns of the new matrices
# are sorted in alphabetical order, by their column name.  If no name appears in
# all of the columns, then NULL is returned.

columnIntersection <- function(L){
  n <- length(L)
  cnames <- colnames(L[[1]])

  for (i in 1:n) {
    cnames <- intersect(cnames, colnames(L[[i]]))
  }

  if (length(cnames) == 0) {
    return(NULL)
  }

  cnames <- sort(cnames)
  instersectL <- vector("list", n)

  for (i in 1:n) {
    instersectL[[i]] <- L[[i]][, match(cnames, colnames(L[[i]])), drop=FALSE]
  }

  return(instersectL)
}

sampleIntersection <- function(data, sampleRows) {
  if (sampleRows) {
    data <- lapply(data, t)
  }
  
  data <- columnIntersection(data)
  
  if (sampleRows) {
    data <- lapply(data, t)
  }
  
  return(data)
}