rep.row <- function (x, n) {
  matrix(rep(x, each=n), nrow=n)
}


groupwise.discover.test.coverage <- function (events, bg) {
  p.covered <- 1 - exp(colSums(log1p(-bg)))
  x <- sum(colSums(events) > 0)
  ppoisbinom(x, p.covered, lower.tail=FALSE)
}


groupwise.discover.test.exclusivity <- function (events, bg) {
  log.p <- log(bg)
  log.not.p <- log1p(-bg)
  p.exactly.one <- exp(matrixStats::colLogSumExps(rep.row(colSums(log.not.p), nrow(events)) + log.p - log.not.p))
  x <- sum(colSums(events) == 1)
  ppoisbinom(x, p.exactly.one, lower.tail=FALSE)
}


groupwise.discover.test.impurity <- function (events, bg) {
  log.p <- log(bg)
  log.not.p <- log1p(-bg)
  p.exactly.one <- exp(matrixStats::colLogSumExps(rep.row(colSums(log.not.p), nrow(events)) + log.p - log.not.p))
  p.none <- exp(colSums(log.not.p))
  p.more.than.one <- 1 - p.none - p.exactly.one
  x <- sum(colSums(events) > 1)
  ppoisbinom(x, p.more.than.one)
}


#' Perform a groupwise mutual exclusivity test.
#'
#' @param x An object of type \code{discover.matrix}, with rows
#'   corresponding to the genes in the gene set to be tested.
#' @param method The mutual exclusivity statistic to estimate
#'   significance for. Possible values are \code{"impurity"},
#'   \code{"coverage"}, and \code{"exclusivity"}.
#' @return An object of type \code{groupwise.discover.out}.
#' 
#' @export
groupwise.discover.test <- function (x, method=c("impurity", "coverage", "exclusivity")) {
  method <- match.arg(method)

  events <- x$events
  bg <- x$bg
  
  grouptest.func <- switch(method,
                           coverage=groupwise.discover.test.coverage,
                           exclusivity=groupwise.discover.test.exclusivity,
                           impurity=groupwise.discover.test.impurity)
  p <- grouptest.func(events, bg)

  result <- list(p.value=p, method=method)
  class(result) <- "groupwise.discover.out"
  result
}
