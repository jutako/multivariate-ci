##########################################################################################
## Plots
##########################################################################################

#' Plot CB when also whole rows have been removed
#'
#' @param dmat [N,M] numeric, Data matrix
#' @param Z [N,M] logical, Output of some mvci method, mvci$Z
#' @param upmask [N,M] logical, Output of some mvci method, mvci$upmask
#' @param downmask [N,M] logical, Output of some mvci method, mvci$downmask
#' @param type char {'inside','outside'}, 'outside' -> plot rows outside CB in red, 'inside' -> plot rows inside CB in red
#' @param title char, Plot title
#' @param xlab char, Plot x axis label
#' @param ylab char, Plot y axis label
#'
#' @return Base graphics plot to current plotting device
#'
#' @importFrom graphics lines plot points
#'
#' @export
plot_data_ci <- function(dmat, Z, upmask, downmask, type = 'outside',
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
#'
#' @param dmat [N,M] numeric, Data matrix
#' @param upmask [N,M] logical, Output of some mvci method, mvci$upmask
#' @param downmask [N,M] logical, Output of some mvci method, mvci$downmask
#' @param timevec [1,M] numeric, Plot positions of the columns of dmat
#' @param title char, Plot title
#' @param xlab char, Plot x axis label
#' @param ylab char, Plot y axis label
#' @param ylim [1,2] numeric, y axis limits
#' @param draw.cb.points [1,1] logical, Draw points that define CB boundary
#' @param samp.frac [1,1] numeric, Sampling fraction from range [0,1]
#' @param yintercept [1,1] numeric, Horizontal line y axis intercept, default: no horizontal line
#'
#' @return Base graphics plot to current plotting device
#'
#' @importFrom graphics lines plot points
#'
#' @export
plot_data_ci2 <- function(dmat, upmask, downmask,
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
#'
#' @param dmat [N,M] numeric, Data matrix
#' @param cb list, The output of some mvci method
#' @param timevec [1,M] numeric, Plot positions of the columns of dmat
#' @param title char, Plot title
#' @param xlab char, Plot x axis label
#' @param ylab char, Plot y axis label
#' @param ylim [1,2] numeric, y axis limits
#' @param xlim [1,2] numeric, x axis limits
#' @param draw.cb.points [1,1] logical, Draw points that define CB boundary
#' @param highlight.rows.out [1,1] logical, Should curves that exit CB be highlighted?
#' @param samp.frac [1,1] numeric, Sampling fraction from range [0,1]
#' @param yintercept [1,1] numeric, Horizontal line y axis intercept, default: no horizontal line
#' @param xaxt char, Passed on to base::plot()
#'
#' @return Base graphics plot to current plotting device
#'
#' @importFrom graphics lines plot points
#'
#' @export
plot_data_cb <- function(dmat, cb,
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
