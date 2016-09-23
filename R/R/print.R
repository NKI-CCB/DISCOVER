#' @export
print.pairwise.discover.out <- function (x, fdr.threshold=0.01) {
  type <- switch(x$alternative, less="mutual exclusivity", greater="co-occurrence")
  cat("Pairwise DISCOVER", type, "test\n")
  cat("alternative hypothesis: observed overlap is", x$alternative, "than expected by chance\n")

  cat("\n")

  cat("number of pairs tested:", sum(!is.na(x$p.values)))
  cat("\n")

  cat("proportion of true null hypotheses:", x$pi0)
  cat("\n")

  cat("number of significant pairs at a maximum FDR of", fdr.threshold, ":", sum(x$pi0 * x$q.values < 0.01, na.rm=TRUE))
  cat("\n")
}


#' @export
print.groupwise.discover.out <- function (x) {
  cat("Groupwise DISCOVER mutual exclusivity test\n")
  cat("\n")
  cat("method:", x$method)
  cat("\n")
  cat("P =", x$p.value)
  cat("\n")
}
