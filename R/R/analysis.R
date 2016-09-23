#' Perform many pairwise mutual exclusivity or co-occurrence tests.
#'
#' @param x A \code{discover.matrix} object.
#' @param g An optional grouping vector for the rows of \code{events}. Pairs of rows within the same group are not tested.
#' @param alternative Either \code{"less"} for mutual-exclusivity analysis, or \code{"greater"} for co-occurrence analysis.
#' @param correct If \code{TRUE}, multiple testing correction is performed.
#' @return An object of type \code{pairwise.discover.out}.
#'
#' @useDynLib discover
#' @export
pairwise.discover.test <- function (x, g=NULL, alternative=c("less", "greater"), correct=TRUE) {
  alternative <- match.arg(alternative)

  events <- x$events
  bg <- x$bg
  
  if (is.null(g)) {
    result <- .Fortran("mutex",
                       as.integer(nrow(events)), as.integer(ncol(events)),
                       as.integer(events), as.double(bg), as.integer(alternative == "less"),
                       p=double(nrow(events) * (nrow(events) - 1) / 2),
                       q=double(nrow(events) * (nrow(events) - 1) / 2),
                       pi0=double(1))
    p <- matrix(NA, nrow(events), nrow(events))
    p[lower.tri(p)] <- result$p

    q <- matrix(NA, nrow(events), nrow(events))
    q[lower.tri(q)] <- result$q
  }
  else {
    i <- order(g)
    block.sizes <- table(g)
    result <- .Fortran("analyseBlockStructure",
                       as.integer(nrow(events)), as.integer(ncol(events)),
                       as.integer(events[i, ]), as.double(bg[i, ]), as.integer(alternative == "less"),
                       as.integer(length(block.sizes)), as.integer(block.sizes),
                       p=double(nrow(events)**2), q=double(nrow(events)**2), pi0=double(1))

    j <- order(i)
    p <- matrix(result$p, nrow=nrow(events))[j, j]
    p[is.nan(p)] <- NA

    q <- matrix(result$q, nrow=nrow(events))[j, j]
    q[is.nan(q)] <- NA
  }

  rownames(p) <- rownames(events)
  colnames(p) <- rownames(events)

  rownames(q) <- rownames(events)
  colnames(q) <- rownames(events)
  
  result <- list(
    p.values=p,
    q.values=q,
    pi0=result$pi0,
    alternative=alternative)
  class(result) <- "pairwise.discover.out"
  result
}


#' @export
as.data.frame.pairwise.discover.out <- function (result, q.threshold=0.01) {
  sig.pairs <- which(result$pi0 * result$q.values < q.threshold, arr.ind=TRUE)
  df <- data.frame(
    gene1=rownames(result$q.values)[sig.pairs[, "row"]],
    gene2=colnames(result$q.values)[sig.pairs[, "col"]],
    p.value=result$p.values[sig.pairs],
    q.value=result$q.values[sig.pairs],
    stringsAsFactors=FALSE)
  df[order(df$p.value), ]
}
