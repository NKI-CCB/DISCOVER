#' Estimate a background matrix.
#'
#' @param events A binary event matrix
#' @param strata A stratification vector, if stratified background estimation is needed
#' @return The estimated background rates
#'
#' @examples
#' x <- matrix(as.integer(runif(100) < 0.3), 10, 10)
#' bg <- estimateBackground(x)
#'
#' @useDynLib discover
estimateBackground <- function (events, strata=NULL) {
  if (is.null(strata)) {
    row.sums <- rowSums(events)
    col.sums <- colSums(events)

    row.values <- unique(row.sums)
    row.inverse <- match(row.sums, row.values)
    row.weights <- table(row.inverse)
    col.values <- unique(col.sums)
    col.inverse <- match(col.sums, col.values)
    col.weights <- table(col.inverse)

    mu <- .Fortran("estimateBackground",
                   as.integer(length(row.values)), as.integer(row.values), as.integer(row.weights),
                   as.integer(length(col.values)), as.integer(col.values), as.integer(col.weights),
                   mu=double(length(row.values) + length(col.values)), PACKAGE="discover")$mu

    num.rows <- length(row.values)
    eA <- exp(mu[1:num.rows] / row.weights) %*% t(exp(mu[(num.rows+1):length(mu)] / col.weights))
    bg <- 1.0 / (eA + 1)

    bg <- bg[row.inverse, col.inverse]
  }
  else {
    bg <- matrix(nrow=nrow(events), ncol=ncol(events))
    for (s in unique(strata)) {
      i <- which(strata == s)
      bg[, i] <- estimateBackground(events[, i])
    }
  }

  dimnames(bg) <- dimnames(events)
  bg
}
