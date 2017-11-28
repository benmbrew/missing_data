.dominateset <- function (xx, KK = 20) {
  zero <- function(x) {
    s = sort(x, decreasing=T, index.return = TRUE)
    x[x < s$x[KK]] = 0
    return(x)
  }
  normalize <- function(X) X/rowSums(X)
  A = matrix(0, nrow(xx), ncol(xx))
  for (i in 1:nrow(xx)) {
    A[i, ] = zero(xx[i, ])
  }
  return(normalize(A))
}