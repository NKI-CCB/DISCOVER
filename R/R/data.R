#' Class to store a binary alteration matrix and the corresponding
#' alteration probability estimates.
#'
#' @param events A binary alteration matrix.
#' @param bg Matrix of alteration probabilities for \code{events}.
#'   Mainly for internal use only. Use at your own risk.
#' @param strata To perform a stratified DISCOVER test, this array
#'   should contain the strata. The length of this array must match the
#'   number of columns of \code{events}.
#' @return An object of type \code{discover.matrix}.
#' 
#' @export
discover.matrix <- function (events, bg=NULL, strata=NULL) {
  events <- as.matrix(events)
  
  if (is.null(bg))
    bg <- estimateBackground(events, strata)
  else if (!is.null(strata))
    warning("The value of 'strata' is ignored if 'bg' is provided.")

  result <- list(events=events, bg=bg)
  class(result) <- "discover.matrix"
  check.consistency(result)
  result
}


check.consistency <- function (x) {
  stopifnot(all(dim(x$data) == dim(x$bg)))
  stopifnot(all(rownames(x$data) == rownames(x$bg)))
  stopifnot(all(colnames(x$data) == rownames(x$bg)))
}


#' @export
`[.discover.matrix` <- function (x, ...) {
  discover.matrix(events=x$events[...], bg=x$bg[...])
}


#' @export
dimnames.discover.matrix <- function (x) {
  dimnames(x$events)
}


#' @export
`dimnames<-.discover.matrix` <- function (x, value) {
  dimnames(x$events) <- value
  dimnames(x$bg) <- value
  check.consistency(x)
  x
}


#' @export
dim.discover.matrix <- function (x) {
  dim(x$events)
}


#' @export
print.discover.matrix <- function (x) {
  cat("DISCOVER data:\n")
  print(x$events)
}


#' @export
#' @method rbind discover.matrix
rbind.discover.matrix <- function (...) {
  matrices <- list(...)
  stopifnot(all(sapply(matrices, function (x) is(x, "discover.matrix"))))
  discover.matrix(events=do.call(rbind, lapply(matrices, function (x) x$events)),
                bg=do.call(rbind, lapply(matrices, function (x) x$bg)))
}


#' @export
#' @method cbind discover.matrix
cbind.discover.matrix <- function (...) {
  matrices <- list(...)
  stopifnot(all(sapply(matrices, function (x) is(x, "discover.matrix"))))
  discover.matrix(events=do.call(cbind, lapply(matrices, function (x) x$events)),
                bg=do.call(cbind, lapply(matrices, function (x) x$bg)))
}
