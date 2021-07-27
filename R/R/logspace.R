#' Compute the logarithm of the sum of exponentials for each column in a matrix.
#'
#' @param x a real-valued matrix.
#' @return a vector with elements corresponding to the log sums of exponentials of the columns of \code{x}.
#'
#' @useDynLib discover
colLogSumExps <- function (x) {
  .Fortran("colLogSumExps", as.double(x), as.integer(nrow(x)), as.integer(ncol(x)), result=double(ncol(x)), NAOK=TRUE)$result
}
