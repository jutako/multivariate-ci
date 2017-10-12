##########################################################################################
## Plots
##########################################################################################

#' Plot CB when also whole rows have been removed
plot.data.ci <- function(dmat, Z, upmask, downmask, type = 'outside',
                         title = 'CI plot',
                         xlab = 'time / variables',
                         ylab = 'value'){

  plot( NA, ylab = ylab, xlab = xlab, main = title,
        xlim=c(1,ncol(dmat)), ylim=c(min(dmat), max(dmat)) )
  for ( i in 1:nrow(dmat) ) {
    lines( dmat[i,], lwd=0.75, col='cyan' )
  }

  row.idx <- which(getRows(Z, type))
  for ( i in row.idx ) {
    lines( dmat[i,], lwd=0.75, col='red' )
  }

  ciup <- dmat[upmask]
  points(1:ncol(dmat), ciup, lwd=1.5, col='blue')
  lines(1:ncol(dmat), ciup, lwd=1.5, col='blue')

  cidown <- dmat[downmask]
  points(1:ncol(dmat), cidown, lwd=1.5, col='blue')
  lines(1:ncol(dmat), cidown, lwd=1.5, col='blue')
}


#' Plot CB if only points have been removed
plot.data.ci2 <- function(dmat, upmask, downmask,
                          timevec = 1:ncol(dmat),
                          title = 'CI plot',
                          xlab = 'time / variables',
                          ylab = 'value',
                          ylim = c(min(dmat), max(dmat)),
                          draw.cb.points = T,
                          samp.frac = 1,
                          yintercept = NULL){

  plot(timevec, colMeans(dmat), lwd=0.75, col='white', ylab = ylab, xlab = xlab, main = title,
       xlim = c(timevec[1], timevec[length(timevec)]), ylim = ylim )

  idx <- sample(1:nrow(dmat), floor(samp.frac * nrow(dmat)), replace = F)
  for ( i in idx ) {
    lines(timevec, dmat[i,], lwd=0.75, col='cyan' )
  }

  ciup <- dmat[upmask]
  lines(timevec, ciup, lwd=1.5, col='blue')

  cidown <- dmat[downmask]
  lines(timevec, cidown, lwd=1.5, col='blue')

  if (draw.cb.points){
    points(timevec, ciup, lwd=1.5, col='blue')
    points(timevec, cidown, lwd=1.5, col='blue')
  }

  if (!is.null(yintercept)){
    lines(timevec, rep(yintercept, length(timevec)), lwd=1.5, col='red')
  }
}


#' Plot CB
plot.data.cb <- function(dmat, cb,
                         timevec = 1:ncol(dmat),
                         title = 'CB plot',
                         xlab = 'time / variables',
                         ylab = 'value',
                         ylim = c(min(dmat), max(dmat)),
                         xlim = c(timevec[1], timevec[length(timevec)]),
                         draw.cb.points = T,
                         highlight.rows.out = T,
                         samp.frac = 1,
                         yintercept = NULL,
                         xaxt = 's'){
  N <- nrow(dmat)

  ciup <- dmat[cb$upmask]
  cidown <- dmat[cb$downmask]
  dmat.out <- dmat[setdiff(1:N, cb$row.inc.idx),]
  dmat <- dmat[cb$row.inc.idx, ]

  # rows inside bands
  plot(timevec, colMeans(dmat), lwd=0.75, col='white', ylab = ylab, xlab = xlab, main = title,
       xlim = xlim, ylim = ylim, xaxt = xaxt)

  idx <- sample(1:nrow(dmat), floor(samp.frac * nrow(dmat)), replace = F)
  for ( i in idx ) {
    lines(timevec, dmat[i,], lwd=0.75, col='cyan' )
  }

  # points outside bands
  if (nrow(dmat.out) > 0){
    if (highlight.rows.out){
      for ( i in 1:nrow(dmat.out) ) {
        lines(timevec, dmat.out[i,], lwd=0.75, col='orange')
      }
    } else {
      for ( i in 1:nrow(dmat.out) ) {
        lines(timevec, dmat.out[i,], lwd=0.75, col='cyan')
      }
    }
  }


  # confidence bands
  lines(timevec, ciup, lwd=2.5, col='blue')
  lines(timevec, cidown, lwd=2.5, col='blue')

  if (draw.cb.points){
    points(timevec, ciup, lwd=1.5, col='blue')
    points(timevec, cidown, lwd=1.5, col='blue')
  }

  # yintercept
  if (!is.null(yintercept)){
    lines(timevec, rep(yintercept, length(timevec)), lwd=1.5, col='black')
  }

}



plot.stockd <- function(cb){

  m <- matrix(c(1,1,1,2), ncol = 1)
  layout(m)
  #layout.show()

  titlestr <- sprintf('stockdata, K=%d, L=%d', cb$K, cb$L)
  plot.data.cb(cb$data, cb,
               timevec = cb$time,
               title = titlestr,
               draw.cb.points = F,
               ylim = c(-200, 800),
               yintercept = 0)

  cs <- colSums(!cb$Z)
  barplot(cs-cb$K)

}


