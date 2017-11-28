standardNormalization <- function (x) 
{
  x = as.matrix(x)
  mean = apply(x, 2, mean, na.rm=TRUE)
  sd = apply(x, 2, sd, na.rm=TRUE)
  sd[sd == 0] = 1
  xNorm = t((t(x) - mean)/sd)
  return(xNorm)
}