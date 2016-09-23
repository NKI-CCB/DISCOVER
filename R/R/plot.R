sort.by.freq <- function (events) {
  gene.order <- order(rowSums(events), decreasing=TRUE)
  sample.order <- rev(do.call(order, as.data.frame(t(events[gene.order, ]))))
  events[gene.order, sample.order]
}


#' @export
plot.discover.matrix <- function(x, sort=TRUE) {
  if(sort)
    events <- sort.by.freq(x$events)
  else
    events <- x$events

  ngenes <- nrow(events)
  nsamples <- ncol(events)
  coverage <- sum(colSums(events) > 0)

  numOfOncos <- ngenes * nsamples
  oncoCords <- matrix(rep(0, numOfOncos * 5), nrow=numOfOncos)
  colnames(oncoCords) <- c("xleft", "ybottom", "xright", "ytop", "altered")

  xpadding <- .01
  ypadding <- .01
  cnt <- 1
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      xleft <- j - 1 + xpadding
      ybottom <- ((ngenes - i + 1) - 1) + ypadding
      xright <- j - xpadding
      ytop <- (ngenes - i + 1) - ypadding
      altered <- events[i, j]
      
      oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
      cnt <- cnt + 1
    }
  }
  
  colors <- rep("lightgray", cnt)
  colors[which(oncoCords[, "altered"] == 1)] <- "black"
  plot(c(0, nsamples), c(0, ngenes), type="n", main=sprintf("Gene set altered in %.2f%%: %d of %d cases", coverage/nsamples*100, coverage, nsamples), xlab="Samples", ylab="", yaxt="n")
  rect(oncoCords[, "xleft"], oncoCords[, "ybottom"], oncoCords[, "xright"], oncoCords[, "ytop"], col=colors, border=NA)
  axis(2, at=(ngenes:1) - 0.5, labels=rownames(events), las=2)
}
