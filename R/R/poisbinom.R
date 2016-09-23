#' Distribution function for the Poisson-Binomial distribution with success probabilities 'probs'.
#'
#' @param q Quantile.
#' @param probs Vector of success probabilities.
#' @param lower.tail If \code{TRUE}, probabilities are P[X <= x], otherwise, P[X > x].
#' @param log.p If \code{TRUE}, probabilities p are given as log(p).
#' @return The cumulative distribution function.
#'
#' @useDynLib discover
ppoisbinom <- function (q, probs, lower.tail=TRUE, log.p=FALSE) {
  if (!lower.tail)
    q <- q - 1
  
  p <- .Fortran("ppoisbinom", as.integer(q), as.integer(length(probs)), as.double(probs), p=double(1))$p

  if (!lower.tail)
    p <- 1 - p

  if (log.p)
    log(p)
  else
    p
}
