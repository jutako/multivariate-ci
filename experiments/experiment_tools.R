
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



#' Experiment 1: Compute Pr[V(x|x_u,x_l)<=L] for various value of L, N_train and B
#'
run.exp1 <- function(cbfun, X_train, X_test, N.train.arr, M.arr, K.prc, L.arr, B){

  N <- nrow(X_train)
  M <- ncol(X_train)

  #   nres <- length(N.train.arr)* length(M.arr)* length(L.arr) * B
  #   resdf <- data.frame(Ntr = vector('integer', nres),
  #                       M = vector('integer', nres),
  #                       K = vector('integer', nres),
  #                       L = vector('integer', nres),
  #                       B = vector('integer', nres),
  #                       prc.in.cbf = vector('numeric', nres),
  #                       prc.in.env = vector('numeric', nres))

  require('foreach')
  options(cores = 3)

  res <- foreach (l = L.arr, .combine = rbind) %dopar% {
    #l <- L.arr[lidx]
    ind <- 1
    nres <- length(N.train.arr)* length(M.arr) * B
    resdf <- data.frame(Ntr = vector('integer', nres),
                        M = vector('integer', nres),
                        K = vector('integer', nres),
                        L = vector('integer', nres),
                        B = vector('integer', nres),
                        prc.in.cbf = vector('numeric', nres),
                        prc.in.env = vector('numeric', nres))
    for (ntr in N.train.arr){
      for (m in M.arr){
        #cat(sprintf('L:%d, ntr: %d, m: %d - bootstrapping ... \n', l, ntr, m))

        k <- floor(K.prc * ntr)
        cbtmpf <- function(X){cbfun(X, K = k, L = l)}
        cidx <- 1:m

        # bootstrap
        for (b in 1:B){
          ridx <- sample(1:N, ntr, replace = F)
          resdf[ind,] <- data.frame(Ntr = ntr, M = m, K = k, L = l, B = b,
                                    prc.in.cbf = prc.inside.cb(cbtmpf, X_train[ridx, cidx], X_test[,cidx], l),
                                    prc.in.env = prc.inside.cb(cb.envelope, X_train[ridx, cidx], X_test[,cidx], l) )
          ind <- ind + 1
        }

        #         bres <- foreach (b = 1:B, .combine = rbind) %do% {
        #           ridx <- sample(1:N, ntr, replace = F)
        #
        #           data.frame(Ntr = ntr, M = m, K = k, L = l, B = b,
        #                      prc.in.cbf = prc.inside.cb(cbtmpf, X_train[ridx, cidx], X_test[,cidx], l),
        #                      prc.in.env = prc.inside.cb(cb.envelope, X_train[ridx, cidx], X_test[,cidx], l) )
        #         }
        #         resdf[ind:(ind+B-1),] <- bres
        #         ind = ind + B
      } #m
    } #ntr

    resdf
  }#l

  res
}


cb_area_test_loader <- function(path,
                                ks = c( 20, 40, 80, 120 ),
                                ls = ls <- c( 0, 1, 2, 5, 10, 20, 50, 100, 600)) {

  foo <- c()
  for ( k in ks ) {
    for ( l in ls ) {
      fname <- sprintf( '%s/stockdata_rawzero_K%d_L%d.rds', path, k, l )
      cb <- readRDS( fname )
      foo <- rbind( foo, c( k, l, cb_area( cb ) ) )
    }
  }
  foo <- as.data.frame( foo )
  names(foo) <- c('K', 'L', 'A')
  foo
}



## call cb_area_test_loader first, give return value as argument to this fnc
cb_area_test <- function( areas ) {
  ks <- unique( areas$K )

  quartz( height=7, width=14 )
  par( mfrow=c(1,2) )

  plot( NA, xlim=c(min(areas$L), max(areas$L)), ylim=c(min(areas$A),max(areas$A)),
        xlab='L', ylab='Area', main='Area as function of L for K in (20,40,80,100)' )
  colors <- c('red', 'green', 'blue', 'magenta', 'cyan')
  for ( i in 1:length(ks) ) {
    k <- ks[i]
    d <- subset( areas, K==k )
    lines( d$L, d$A, col=colors[i], type='b' )
  }

  plot( NA, xlim=c(min(areas$L), max(areas$L)), ylim=c(0.85,1),
        xlab='L', ylab='area/areaL0',
        main='Reduction in area wrt area of L=0')
  for( i in 1:length(ks) ) {
    k <- ks[i]
    d <- subset( areas, K==k )
    d$A <- d$A/d$A[1]
    lines( d$L, d$A, col=colors[i], type='b' )
  }
}
