
# Implied packages:
require(stats)
require(signal)
require(zoo)


#-----------------------------------------------------------------------------------------
## dense_area_motivation.R

#' Make a dataset to showcase the benefit of having L>0
make_dam_data_2d <- function(N0, out.prc = 0.2,
                             mu0 = 0, sd0 = 1,
                             mu1 = 0, sd1 = 6){
  ## Bulk of data
  x1 <- rnorm(N0, mean = mu0, sd = sd0)
  x2 <- rnorm(N0, mean = mu0, sd = sd0)
  d <- matrix(c(x1,x2), byrow = F, ncol = 2)

  ## Outliers
  N1 <- round(out.prc * N0)
  outl <- rnorm(N1, mean = mu1, sd = sd1)
  idx <- sample(1:N0, N1, replace = F)
  idx1 <- idx[1:floor(N1/2)]
  idx2 <- idx[(floor(N1/2)+1):N1]

  # Add to data
  d[idx1, 1] <- d[idx1, 1] + outl[1:floor(N1/2)]
  d[idx2, 2] <- d[idx2, 2] + outl[(floor(N1/2)+1):N1]

  daf <- as.data.frame(d)
  daf$is.outlier <- F
  daf$is.outlier[idx] <- T

  list(d = d, df = daf)
}


#' Make a dataset to showcase the benefit of having L>0
make_dam_data <- function(N0, M,
                          N.outliers = 0.2*N0,
                          outlier.dim = NULL,
                          mu0 = 0, sd0 = 1,
                          mu1 = 0, sd1 = 6){
  #   M <- 4
  #   N0 <- 100
  #   out.prc = 0.2
  #   mu0 = 0
  #   sd0 = 1
  #   mu1 = 0
  #   sd1 = 6

  ## Bulk of data
  xl <- list()
  for (m in 1:M){
    xl[[m]] <- rnorm(N0, mean = mu0, sd = sd0)
  }
  d <- matrix(unlist(xl), byrow = F, ncol = M)

  ## Outliers
  N1 <- round(N.outliers)
  outl <- rnorm(N1, mean = mu1, sd = sd1)
  idx <- sample(1:N0, N1, replace = T)
  #with replacement to allow several outiers per observation

  # Add to data
  if (is.null(outlier.dim)){
    m.arr <- 1:M
  } else {
    m.arr <- outlier.dim
  }

  for (i in 1:N1){
    if (length(m.arr)>1){
      cur.m <- sample(m.arr, 1)
    } else {
      cur.m <- m.arr
    }

    d[idx[i], cur.m] <- d[idx[i], cur.m] + outl[i]
  }

  daf <- as.data.frame(d)
  daf$is.outlier <- F
  daf$is.outlier[idx] <- T

  list(d = d, df = daf)
}


#' Make smooth "random normal" data vectors
rnorm.smooth <- function(N, mean = 0, sd = 1){
  d <- stats::smooth.spline(rnorm(N+20, mean = mean, sd = sd), spar = 0.6)$y
  d[11:(length(d)-10)]
}


#' Make smooth looking "bump"
#' plot(make.bump.smooth(100), type = 'l')
make.bump.smooth <- function(N, k=floor(0.3*N)){
  set.seed(1) #1 gives a bump
  d <- rnorm(N)
  d <- zoo::rollmean(d, k)
  d <- d * signal::hamming(N-k+1)
  d <- signal::resample(d, p = N, q = length(d))
  d <- stats::smooth.spline(d)$y
  d
}


#' Make a toy dataset
#' matplot(t(make.toy.data1(100, 20, 100)), type = 'l')
make.toy.data1 <- function(Nbase, Nbump, M){

  dmat <- matrix(NA, nrow = Nbase + Nbump + 1, M)
  bump <- make.bump.smooth(M)
  for (i in 1:Nbump){
    dmat[i,] <- bump + 0.1*rnorm.smooth(M) + rnorm(1, sd = 0.05)
  }
  for (i in (Nbump+1):(Nbump+Nbase)){
    dmat[i,] <-  0.1*rnorm.smooth(M)
  }

  dmat[Nbump + Nbase + 1, ] <- 0.01*rnorm.smooth(M)
  dmat[Nbump + Nbase + 1, 10] <- 0.25

  dmat
}


#' Make normal pdf looking "bump"
#' plot(make.bump.normal(10), type = 'l')
make.bump.normal <- function(N){
  s <- dnorm(-4:4, mean = 0, sd = 1.5)
  to.pad.l <- floor((N - length(s))/2)
  to.pad.r <- N - length(s) - to.pad.l
  c(rep(0, to.pad.l), s, rep(0, to.pad.r))
}


#' Make a toy dataset
#' matplot(t(make.toy.data2(50, 50, 40)), type = 'l')
make.toy.data2 <- function(Nbase, Nbump, M, jitter = floor(M/7)){
  #jitter = floor(M/4)

  dmat <- matrix(NA, nrow = Nbase + Nbump + 2, M)
  bump <- make.bump.normal(M)

  # make base
  for (i in 1:Nbase){
    dmat[i,] <-  0.1*rnorm.smooth(M)
  }

  # make bumps
  for (i in (Nbase + 1):(Nbase + Nbump)){
    dmat[i,] <- shifter(bump, sample(-jitter:jitter,1)) + 0.1*rnorm.smooth(M) #+ rnorm(1, sd = 0.05)
  }

  # global outlier
  bump2 <- make.bump.normal(10)
  bump2 <- signal::resample(bump2, p = M , q = length(bump2) )
  bump2 <- 1.1 * c( rev(head(bump2, floor(M/2))), rev(tail(bump2, floor(M/2))) )
  dmat[Nbase + Nbump + 1, ] <- stats::smooth.spline(bump2)$y + 0.025

  # suppressed normal
  dmat[Nbase + Nbump + 2, ] <- 0.02*rnorm.smooth(M)

  class.vec <- rep("", nrow(dmat))
  class.vec[1:Nbase] <- 'normal'
  class.vec[(Nbase + 1) : (Nbase + Nbump)] <- 'local.outlier'
  class.vec[Nbase + Nbump + 1] <- 'global.outlier'
  class.vec[Nbase + Nbump + 2] <- 'supp.normal'

  attr(dmat, 'groups') <- class.vec
  dmat
}


#' Read confidence band from result list
cbr_extract_cb <- function(X, cb){
  matrix(c(X[cb$downmask], X[cb$upmask]), nrow = 2, byrow = T)
}


#' Compute percentage of obsevations that fall within confidence bands using some L criterion
cbr_prc_inside <- function(Xtest, cbm, L){
  total.contained(Xtest, cbm[1,], cbm[2,], L = L)
}


#' Circularly shift vector elemens
#' To e.g. change the position of a peak within a series
#' positive n shifts things to the left, negative to the right
#'
#' From: http://stackoverflow.com/questions/30542128/circular-shifting-arrays-in-r-by-distance-n
shifter <- function(x, n) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}


#' Transform a matrix into a ggplottable data.frame
make.ggplot.df <- function(dmat){
  pd <- as.data.frame(dmat)
  if ( !is.null(dimnames(dmat)) ){
    pd$row <- dimnames(dmat)[[1]]
  } else {
    pd$row <- 1:nrow(pd)
  }
  pd <- reshape2::melt(pd, id.vars = c('row'))
  pd
}

#' Plot data and CB: one of many visualizations
#'
toyplot <- function(dmat, cb_naive, cbm0, cbml){

  # pd <- make.ggplot.df(dmat)
  # pd$id <- 'data'
  # pdn <- make.ggplot.df(matrix(c(cb_naive$cb.low, cb_naive$cb.high), nrow = 2, byrow = T))
  # pdn$id <- 'naive'
  # pd0 <- make.ggplot.df(cbm0)
  # pd0$id <- 'L=0'
  # pdl <- make.ggplot.df(cbml)
  # pdl$id <- 'L>0'
  # pd2 <- rbind(pdn, pd0, pdl)

  # Does not work correctly:
  # p <- ggplot(data = pd, mapping = aes(x = variable, y = value,
  #                                      group = row, color = id))
  # p <- p + geom_line(alpha = 0.5)
  # p <- p + geom_line(data = pd2, size = 1)
  # p <- p + theme_light()
  # p

  pd <- make.ggplot.df(dmat)
  pdn <- make.ggplot.df(matrix(c(cb_naive$cb.low, cb_naive$cb.high), nrow = 2, byrow = T))
  pd0 <- make.ggplot.df(cbm0)
  pdl <- make.ggplot.df(cbml)

  p <- ggplot(data = pd, mapping = aes(x = variable, y = value, group = row))
  p <- p + geom_line(alpha = 0.5)
  p <- p + geom_line(data = pdn, color = 'red', size = 1)
  p <- p + geom_line(data = pd0, color = 'blue', size = 1)
  p <- p + geom_line(data = pdl, color = 'green', size = 1)
  p <- p + theme_light()
  p

}

