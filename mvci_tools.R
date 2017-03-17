
# Implied packages:
# require(foreach) #run.exp1()

##########################################################################################
## Confidence regions
##########################################################################################

#' Find shortest interval quantiles
#'
#' @description
#' Excludes k extreme values from x such that the remaining contiguous range has the
#' shortest width.
#'
#' Note: If the distribution of x is uniform, e.g. x = 1:N, or otherwise there are
#' sevaral shortest intervals then the first is returned, in this example 1:(N-k).
#'
#' @param x [1,N] numeric, A 1D data vector
#' @param k [1,1] int, Number of values to exclude
#'
#' @return [1,2] numeric, The endpoints [min, max] of the quantile
min_width_quantile <- function(x, k) {
  n <- length(x)
  x <- sort(x)
  i <- which.min(x[(n-k):n]-x[1:(k+1)])
  c(x[i],x[n+i-k-1])
}


#' Find central quantiles (densest region of data)
#'
#' @description
#' Symmetrically excludes k extreme values from x, k/2 from each end.
#' If k is odd, then one less value is excluded from the left tail.
#'
#' @param x [1,N] numeric, A 1D data vector
#' @param k [1,1] int, Number of values to exclude
#'
#' @return [1,2] numeric, The endpoints [min, max] of the quantile
central_quantile <- function(x, k){
  n <- length(x)
  x <- sort(x)
  idx <- floor(k/2)
  if (k %% 2 == 0){ #even
    c(x[idx+1],x[n-(k/2)])
  } else { #odd
    i <- floor(k/2)
    c(x[idx+1],x[n-i-1])
  }
}


#' Quantile based naive confidence areas
#'
#' @description
#' No multiplicity corrections.
#'
#' @param dmat [N,M] numeric matrix, N samples, M variables.
#' @param k [1,1] integer, Number of observations to exclude
#' @param quantile.fun function, A function to compute the 1D quantiles with,
#'          defaults to min_width_quantile()
#'          Should accept two arguments: a data vector x and number of observations to
#'          exluce k. See min_width_quantile() and central_quantile() for examples.
#'
#' @return A list with elements:
#' \item{$cb:}{[2,M] numeric array containing the lower and upper confidence band}
#' \item{$k:}{[1,1] integer, the input k for reference}
#' \item{$quantile.fun:}{function, the input quantile.fun for reference}
#' \item{$wall_time_taken:}{[1,1] double, Time taken to find the result, in seconds}
cb_quantile <- function(dmat, k, quantile.fun = min_width_quantile){
  start_time <- Sys.time()

  cb <- apply(dmat, 2, quantile.fun, k = k)

  list( cb = cb,
        k = k,
        quantile.fun = quantile.fun,
        wall_time_taken = as.numeric(difftime(Sys.time(), start_time, units = "sec")) )
}


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



##########################################################################################
## Cross-validation based control procedures
##########################################################################################

#' Generic B -fold crossvalidation curve
#'
#' Input:
#'   cbfun         Function to compute the confidence bands with, interface: cbmat = cbfun(data, k)
#'   cb.contains.fun     Function to test if vector x belongs to a given CB
#'   data          (N,M) numeric matrix, N samples, M dimensions/variables. The samples define the
#'                 empirical cdfs of the m variables.
#'   B             (1,1) integer, How many cross validation folds to use. l = [2,...,N].
#'   k.max  (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#'                 of the test fold observations are outside confidence band the profile computation
#'                 ends. Can be used to speed up computations when only the beginning of the profile
#'                 curve is of interest.
#'
#' Output:
#'   A list with elements:
#'   profile   (1,N-floor(N/B)) numeric, Crossvalidation based profile of FWER control
#'
fold.profile <- function(cbfun, cb.contains.fun, data, B, k.max=floor(0.5*((B-1)/B)*nrow(data)) ){

  N = nrow(data)
  M = ncol(data)

  # Create folds
  data <- data[sample(1:N),] #make sure data is randomized
  fold.len <- floor(N/B)
  fold.idx <- matrix(c(seq(1,N,fold.len)[1:B], seq(fold.len,N,fold.len)[1:B-1], N), nrow=B, ncol=2)
  # (N,2) matrix of fold index [start, stop]
  #fold.idx[,2]-fold.idx[,1]+1

  all.idx <- 1:N
  B.profiles <- matrix(0, nrow=B, ncol=(N-fold.len))
  for (l in 1:B){
    cat(sprintf("round:%d...\n",l))
    if (fold.len==1){
      # Selecting a single row from matrix results in a numeric vector! Casting back to matrix flips
      # the dimensions, thus the transpose. Who the fuck invented this change of class ...
      data.test <- as.matrix(t(data[fold.idx[l,1]:fold.idx[l,2],]))
    } else {
      data.test <- data[fold.idx[l,1]:fold.idx[l,2],]
    }
    #browser()
    B.profiles[l,] <- k.profile(cbfun, cb.contains.fun,
                                data[all.idx[-seq(fold.idx[l,1],fold.idx[l,2],1)],],
                                data.test,
                                k.max = k.max)$k.profile[1:(N-fold.len)]
  }
  B.profile <- colSums(B.profiles, na.rm=T)/N

  #   plot(B.profile, ylim=c(0,1), type="l")
  #   lines(matrix(c(0,0,245,1),nrow=2,byrow=T))
  #   lines(matrix(c(0,alpha,245,alpha),nrow=2,byrow=T))
  #   browser()
  return(list(profile = B.profile,
              profile.mat = B.profiles,
              profile.k = seq(-1,length(B.profile)-2,1),
              fold.idx = fold.idx))
}


#' Generic test set "k profile" (using whichever CB algorithm)
#'
#' Determines the iterations k at which observations in the test set become extreme
#' i.e. do not belong to the confidence band anymore. The result is an (1,N) vector
#' of integers {0,1,...,nrow(data.test} identifying how many of the test set observations
#' are extreme at a given k=[1,...,N].
#'
#' This funcition is needed in the fold.profile() algorithm.
#'
#' Input:
#'   cbfun         Function to compute the confidence bands with, interface: cbmat = cbfun(data, k)
#'   cb.contains.fun     Function to test if vector x belongs to a given CB, TRUE = inside, F = outside
#'   data.train         (N,M) numeric matrix, N samples, M dimensions/variables. The samples define the
#'                   empirical cdfs of the m variables.
#'   data.test       (N.test,M) numeric matrix, Test dataset
#'   k.max    (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#'                   of the observations in 'data.test' are outside confidence band the profile
#'                   computation ends. Can be used to speed up computations when only the beginning of the profile
#'                   curve is of interest.
#'
#' Output:
#'   A list with elements:
#'   k.profile   (1,N) vector of integers [0,1,...,nrow(data.test)]
k.profile <- function(cbfun, cb.contains.fun, data.train, data.test, k.max=floor(0.5*nrow(data.train))){

#   data.train <- makeX(20, 50)
#   data.test <- makeX(5, 50)
#   cbfun <- function(dmat, k){
#     cb <- findcb_topdown(dmat, K = k, L = 1)
#     matrix(c(dmat[cb$downmask], dmat[cb$upmask]), nrow = 2, byrow = T)
#   }
#   k.max = 10
#   cbfun(data.train, 2)

  M = ncol(data.train)
  N.train = nrow(data.train) #size of "training" set
  N.test = nrow(data.test) #size of test set

  # Initialize data structures
  profile.mat <- matrix(1, nrow=N.test, ncol=N.train+2) #store a [0,1] vector for each observation,
  # reserve the 1st column for the case k=0 i.e. data envelope

  # Check data envelope -> store to position 1
  cbe <- cb.envelope(data.train)
  profile.mat[,1] <- as.integer(apply( data.test, 1,
                                       function(x) !cb.contains.fun( x, cbe$envelope[1,], cbe$envelope[2,]) ))
  #profile.mat[,1] <- as.integer(row.outside.ci(data.test, cbe$envelope[1,], cbe$envelope[2,] ))

  continue <- T
  k <- 2
  while (continue){
    cbmat <- cbfun(data.train, k-2)
    profile.mat[,k] <- as.integer(apply( data.test, 1,
                                         function(x) !cb.contains.fun( x, cbmat[1,], cbmat[2,]) ))
    #profile.mat[,k] <- as.integer(row.outside.ci(data.test, cbmat[1,], cbmat[2,]))

    #     debug plot:
    #     matplot(t(matrix(c(cbmat[1,], cbmat[2,], byrow=T, nrow=2)), type="l")
    #     matlines(data.train[rmng.row.arr,], lty=1, lwd=3, col=rgb(0,0,1,1))

    if ((k > k.max + 2) | (sum(profile.mat[,k]) == nrow(profile.mat)) ){
      # k.max has been reached or all test observations are outside CI
      continue = F
    }
    k <- k + 1
  } #of while

  return(list(k.profile = colSums(profile.mat),
              k.profile.mat = profile.mat))
}



## Minimum width envelope confidence region, FWER controlled, l fold crossvalidation
#
# Using cross validation approach aims at approximating the effective k to be used in the mwe.ci()
# function so that FWER is controlled at level alpha. Improves upon mwe.fwer.ci() which for small N
# outputs instabile results. Crossvalidation smooths out the instability.
#
# Input:
#   data          (N,M) numeric matrix, N samples, M dimensions/variables. The samples define the
#                                       empirical cdfs of the m variables.
#   B             (1,1) integer, How many cross validation folds to use. l = [2,...,N].
#   k.max  (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#                 of the test fold observations are outside confidence band the profile computation
#                 ends. Can be used to speed up computations when only the beginning of the profile
#                 curve is of interest.
#
# Output:
#   A list with elements:
#   profile   (1,N-floor(N/B)) numeric, Crossvalidation based profile of FWER control
mwe.fold.profile <- function(data, B, k.max = floor(0.5*((B-1)/B)*nrow(data)) ){

  N = nrow(data)
  M = ncol(data)
  #browser()
  # Create folds
  data <- data[sample(1:N),] #make sure data is randomized
  fold.len <- floor(N/B)
  fold.idx <- matrix(c(seq(1,N,fold.len)[1:B], seq(fold.len,N,fold.len)[1:B-1], N), nrow=B, ncol=2)
  # (N,2) matrix of fold index [start, stop]
  #fold.idx[,2]-fold.idx[,1]+1

  all.idx <- 1:N
  B.profiles <- matrix(0, nrow=B, ncol=(N-fold.len))
  for (b in 1:B){
    cat(sprintf("round:%d...\n",b))
    if (fold.len==1){
      # Selecting a single row from matrix results in a numeric vector! Casting back to matrix flips
      # the dimensions, thus the transpose. Who the fuck invented this change of class ...
      data.test <- as.matrix(t(data[fold.idx[b,1]:fold.idx[b,2],]))
    } else {
      data.test <- data[fold.idx[b,1]:fold.idx[b,2],]
    }
    #browser()
    B.profiles[b,] <- mwe.k.profile(data[all.idx[-seq(fold.idx[b,1],fold.idx[b,2],1)],],
                                    data.test,
                                    k.max)$k.profile[1:(N-fold.len)] / nrow(data.test)
  }

  # todo: how to average over folds
  B.profile <- apply(B.profiles, 2, mean, na.rm=T) #coverage: % of rows inside conf band
  #B.profile <- 1 - colSums(B.profiles, na.rm=T)/N #coverage: % of rows inside conf band
  #browser()
  #   plot(B.profile, ylim=c(0,1), type="l")
  #   lines(matrix(c(0,0,245,1),nrow=2,byrow=T))
  #   lines(matrix(c(0,alpha,245,alpha),nrow=2,byrow=T))
  #   browser()
  return(list(profile = 1 - B.profile,
              profile.mat = B.profiles,
              profile.k = seq(0,length(B.profile)-1,1),
              fold.idx = fold.idx))
}


## Test set "k profile" using Minimum Width Envelope greedy algorithm
#
# Determines the iterations k at which observations in the test set become extreme
# i.e. do not belong to the confidence band anymore. The result is an (1,N) vector
# of integers {0,1,...,nrow(data.test} identifying how many of the test set observations
# are extreme at a given k=[1,...,N].
#
# This funcition is needed in the mwe.fold() algorithm.
#
# Input:
#   data.rs         (N,M) numeric matrix, N samples, M dimensions/variables. The samples define the
#                   empirical cdfs of the m variables.
#   data.test       (N.test,M) numeric matrix, Test dataset
#   k.max    (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#                   of the observations in 'data.test' are outside confidence band the profile
#                   computation ends. Can be used to speed up computations when only the beginning of the profile
#                   curve is of interest.
#
# Output:
#   A list with elements:
#   k.profile   (1,N) vector of integers [0,1,...,nrow(data.test)]
mwe.k.profile <- function(data.rs, data.test, k.max){

  N.mwe = nrow(data.rs) #size of MWE computation set
  M = ncol(data.rs)
  N.test = nrow(data.test) #size of test set

  # (1) Initialize data structures:
  profile.mat <- matrix(1, nrow=N.test, ncol=N.mwe+1) #store a [0,1] vector for each observation,
  # reserve the 1st column for the case k=0 i.e. data envelope

  # 'data.rank' = [N.mwe,M] matrix of columnwise ranks, for each column separately
  data.rank <- apply(data.rs, 2, function(x){order(order(x))}) #(MNlgN)

  # 'CER' = "current extreme ranks" structure
  # do this in a simpler way (order)
  CER <- list(r.max = matrix(c(rep(N.mwe,M),rep(N.mwe-1,M)),nrow=2, byrow=T),
              n.max = matrix(c(apply(data.rank, 2, function(x){which(x==N.mwe)}),
                               apply(data.rank, 2, function(x){which(x==(N.mwe-1))})), nrow=2, byrow=T),
              r.min = matrix(c(rep(1,M),rep(2,M)),nrow=2, byrow=T),
              n.min = matrix(c(apply(data.rank, 2, function(x){which(x==1)}),
                               apply(data.rank, 2, function(x){which(x==(2))})), nrow=2, byrow=T))
  # Store current envelope (using already extracted indices)
  vec.idx.max = 0:(ncol(data.rs)-1)*nrow(data.rs)+CER$n.max[1,]
  vec.idx.min = 0:(ncol(data.rs)-1)*nrow(data.rs)+CER$n.min[1,]
  CER <- c(CER, list(env.max=data.rs[vec.idx.max], env.min=data.rs[vec.idx.min]))
  profile.mat[,1] <- as.integer(row.outside.ci(data.test, CER$env.min, CER$env.max))

  # O(MN)

  n.extr.rows <- 0
  extreme.obs.inds <- NULL
  # Keep track of which rows have not yet been removed
  rmng.row.arr <- 1:N.mwe
  continue <- T
  k <- 2
  while (continue){

    # (2) Find the optimal row to remove (maximize confidence region shrinkage)
    dU <- rep(0,N.mwe)

    for (m in 1:M){
      # Add gains from extreme points to the corresponding observations
      dU[CER$n.max[1,m]] <- dU[CER$n.max[1,m]] + abs(data.rs[CER$n.max[1,m],m]-data.rs[CER$n.max[2,m],m])
      dU[CER$n.min[1,m]] <- dU[CER$n.min[1,m]] + abs(data.rs[CER$n.min[1,m],m]-data.rs[CER$n.min[2,m],m])
    } # O(M)

    n.opt <- which.max(dU) #O(N.mwe)
    extreme.obs.inds <- c(extreme.obs.inds, n.opt)
    n.extr.rows <- n.extr.rows + 1

    # (3) Cleanup & book keeping for the next round
    rmng.row.arr <- setdiff(rmng.row.arr, n.opt) #O(?)

    for (m in 1:M){
      if (CER$n.min[1,m] == n.opt){
        CER$n.min[1,m] = CER$n.min[2,m]
        CER$r.min[1,m] = CER$r.min[2,m]
        CER$env.min[m] = data.rs[CER$n.min[1,m],m]

        min.df <- min.n.set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
        #Note: rmng.row.arr contains the previous 2nd smallest row -> hence look for 2nd smallest el
        CER$n.min[2,m] = rmng.row.arr[min.df$ind[2]]
        CER$r.min[2,m] = min.df$value[2]

      } else if (CER$n.min[2,m] == n.opt){
        min.df <- min.n.set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
        CER$n.min[2,m] = rmng.row.arr[min.df$ind[2]]
        CER$r.min[2,m] = min.df$value[2]

      } else if (CER$n.max[1,m] == n.opt){
        CER$n.max[1,m] = CER$n.max[2,m]
        CER$r.max[1,m] = CER$r.max[2,m]
        CER$env.max[m] = data.rs[CER$n.max[1,m],m]

        max.df <- max.n.set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
        CER$n.max[2,m] = rmng.row.arr[max.df$ind[2]]
        CER$r.max[2,m] = max.df$value[2]

      } else if (CER$n.max[2,m] == n.opt){
        max.df <- max.n.set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
        CER$n.max[2,m] = rmng.row.arr[max.df$ind[2]]
        CER$r.max[2,m] = max.df$value[2]
      }
    } # O(MN)

    # Check if new test observations have become extreme
    profile.mat[,k] <- as.integer(row.outside.ci(data.test, CER$env.min, CER$env.max))

    #     debug plot:
    #     matplot(t(matrix(c(CER$env.min, CER$env.max), byrow=T, nrow=2)), type="l")
    #     matlines(data.rs[rmng.row.arr,], lty=1, lwd=3, col=rgb(0,0,1,1))

    if (k.max <= k){
      continue = F # the desired alpha level k/N has been reached
    }

    k <- k + 1
  } #of while

  # returns numbers and matches of rows that are _outside_ conf bands at a given k
  return(list(k.profile = colSums(profile.mat),
              k.profile.mat = profile.mat))
}


# An extension of max()
#
# Searches for maximum values up to k:th largest value. See k.min() for details.
max.n.set <- function(x,k){
  if (k <= 0){ stop("k.min: too small k") }
  if (length(x) < k){ stop("k.min: too large k") }
  res <- min.n.set(-x,k)
  res$value <- -res$value
  return(res)
}

# Find n:th largest value of x
# From: http://stackoverflow.com/questions/2453326/fastest-way-to-find-second-third-highest-lowest-value-in-vector-or-column
max.n <- function(x, n){
  len <- length(x)
  if(n > len){
    warning('N greater than length(x).  Setting N=length(x)')
    n <- length(x)
  }
  sort(x, partial = len-n+1)[len-n+1]
}

# An extension of min()
#
# Searches for minimum values up to k:th smallest value in O(k*N) time.
#
# Input:
# x   (1,N) numeric, Data vector
# k   (1,1) integer, Up to how far to search
#
# Output:
# res data.frame with k rows and columns:
#     k: 1 for min, 2 for 2nd smallest, etc. ...
#     ind: position within x
#     value: x[ind]
#     so that x[res$ind] = res$value are the k smallest values, starting with the smallest
min.n.set <- function(x,k){
  if (k <= 0){ stop("k.min: too small k") }
  if (length(x) < k){ stop("k.min: too large k") }

  res <- data.frame()
  df <- data.frame(ind=1:length(x), value=x)
  for (i in 1:k){
    pos <- which.min(df$value)
    res <- rbind(res, data.frame(k=i, ind=df$ind[pos], value=df$value[pos]))
    df <- df[-pos,]
    rm(pos)
  }
  return(res)
}

# Find n:th smallest value of x
min.n <- function(x, n){
  -max.n(-x, n)
}


##########################################################################################
## Misc
##########################################################################################


## Observation envelope and related statistics
#
# Inputs:
#   data  [N,M] numeric matrix, N samples, M variables/dimensions.
#
# Outputs:
#   list  envelope    [2,M] matrix, lower (row 1) and upper (row 2) envelopes
#         sum.width   [1,1] numeric, Envelope width (sum of widths over the M dimensions)
#         mean.width  [1,1] numeric, Mean envelope width (mean of widths over the M dimensions)
#         min/max     you get it...
#
cb.envelope <- function(data){
  N <- nrow(data)
  M <- ncol(data)

  order.matrix <- apply(data, 2, order)

  downmask <- matrix(F, nrow = N, ncol = M)
  upmask <- matrix(F, nrow = N, ncol = M)
  for (m in 1:M){
    downmask[order.matrix[1,m] ,m] <- T
    upmask[order.matrix[N,m] ,m] <- T
  }

  envelope <- matrix(c(data[downmask], data[upmask]), nrow=2, byrow=T)
  env.diff = abs(envelope[2,]-envelope[1,])

  list(downmask = downmask,
      upmask = upmask,
      envelope=envelope,
      sum.width=sum(env.diff),
      mean.width=mean(env.diff),
      max=max(env.diff),
      min=min(env.diff))
}


# Find observations that are outside confidence intervals in at least one variable
#
# Inputs:
#   data.train   (n,m) numeric matrix, n samples, m variables. The samples define the empirical cdfs of the m variables.
#             Variables correspond to multiple hypotheses.
#   cb.low    (1,m) numeric vector, lower CI limit for each of the dimensions
#   cb.high   (1,m) numeric vector, upper CI limit for each of the dimensions
row.outside.ci <- function(data.train, cb.low, cb.high){
  N = nrow(data.train)
  data.out.match <- (data.train < repmat(matrix(cb.low, nrow=1, byrow=T), N, 1)) |
    (repmat(matrix(cb.high, nrow=1, byrow=T), N, 1) < data.train)
  return(0 < rowSums(data.out.match))
}


## repmat (as in Matlab)
# http://waxworksmath.com/Authors/G_M/Hastie/Code/Chapter4/repmat.R
repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

#' Count of cb violations in vector x
cb.violations <- function( x, cb.low, cb.high ) {
  sum( sapply( 1:length(cb.high), function(i) cb.low[i] > x[i] || cb.high[i] < x[i] ) )
}

#' Does vector x violate cb with criterion L
cb.contains <- function( x, cb.low, cb.high, L ) {
  cb.violations( x, cb.low, cb.high ) <= L
}

#' Total number of rows in X that do not violate cb
total.contained <- function( X, cb.low, cb.high, L ) {
  sum( apply( X, 1, function(x) cb.contains( x, cb.low, cb.high, L ) ) )
}


#' Compute Pr[V(x|x_u,x_l)<=L] using X_train, K, L and X_test
#'
prc.inside.cb <- function(cbfun, X_train, X_test, L){

# cbfun <- findcb_topdown
# X_train <- makeX(100, 20)
# X_test <- makeX(100, 20)
# K = 5
# L = 1

  cb <- cbfun(X_train)
  in.count <- sum(apply( X_test, 1,
                    function(x) cb.contains( x, X_train[cb$downmask], X_train[cb$upmask], L) ))
  in.count/nrow(X_test)
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


cb_area <- function( cb ) {
  ciup   <- cb$data[cb$upmask]
  cidown <- cb$data[cb$downmask]

  sum( ciup - cidown )
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

