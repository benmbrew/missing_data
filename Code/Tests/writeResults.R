# Append the vector of (argv, method, results) to a new line in file.
# results: vector of results to be saved
# file: path to file which will contain the results
writeResults <- function(results, file) {
  write(results, file, ncolumns=length(results), append=TRUE)
}