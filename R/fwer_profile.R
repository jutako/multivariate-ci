##########################################################################################
## Cross-validation based control procedures
##########################################################################################

#' Generic B -fold crossvalidation curve
#'
#' @param cbfun Function to compute the confidence bands with, interface: cbmat = cbfun(data, k)
#' @param cb.contains.fun Function to test if vector x belongs to a given CB
#' @param data (N,M) numeric matrix, N samples, M dimensions/variables.
#'             The samples define the empirical cdfs of the M variables.
#' @param B (1,1) integer, How many cross validation folds to use. l = [2,...,N].
#' @param k.max (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#'                 of the test fold observations are outside confidence band the profile computation
#'                 ends. Can be used to speed up computations when only the beginning of the profile
#'                 curve is of interest.
#'
#' @return A list with elements:
#' \item{profile}{(1,N-floor(N/B)) numeric, Crossvalidation based profile of FWER control}
#'
#' @export
fold.profile <- function(cbfun, cb.contains.fun, data, B,
                         k.max=floor(0.5*((B-1)/B)*nrow(data)) ){

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
#' @description
#' Determines the iterations k at which observations in the test set become extreme
#' i.e. do not belong to the confidence band anymore. The result is an (1,N) vector
#' of integers {0,1,...,nrow(data.test} identifying how many of the test set observations
#' are extreme at a given k=[1,...,N].
#'
#' This function is needed in the fold.profile() algorithm.
#'
#' @param cbfun Function to compute the confidence bands with, interface: cbmat = cbfun(data, k)
#' @param cb.contains.fun Function to test if vector x belongs to a given CB, TRUE = inside, F = outside
#' @param data.train (N,M) numeric matrix, N samples, M dimensions/variables. The samples define the
#'                   empirical cdfs of the m variables.
#' @param data.test (N.test,M) numeric matrix, Test dataset
#' @param k.max (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#'                   of the observations in 'data.test' are outside confidence band the profile
#'                   computation ends. Can be used to speed up computations when only the beginning of the profile
#'                   curve is of interest.
#'
#' @return A list with elements:
#' \item{k.profile}{(1,N) vector of integers [0,1,...,nrow(data.test)]}
#'
#' @export
k.profile <- function(cbfun, cb.contains.fun, data.train, data.test,
                      k.max=floor(0.5*nrow(data.train))){

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



#' Minimum width envelope confidence region, FWER controlled, l fold crossvalidation
#'
#' @description
#' Using cross validation approach aims at approximating the effective k to be used in the mwe.ci()
#' function so that FWER is controlled at level alpha. Improves upon mwe.fwer.ci() which for small N
#' outputs instabile results. Crossvalidation smooths out the instability.
#'
#' @param data (N,M) numeric matrix, N samples, M dimensions/variables.
#'             The samples define the empirical cdfs of the m variables.
#' @param B (1,1) integer, How many cross validation folds to use. l = [2,...,N].
#' @param k.max (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#'                 of the test fold observations are outside confidence band the profile computation
#'                 ends. Can be used to speed up computations when only the beginning of the profile
#'                 curve is of interest.
#'
#' @return A list with elements:
#' \item{profile}{(1,N-floor(N/B)) numeric, Crossvalidation based profile of FWER control}
#'
#' @export
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
    # profiles according to style a) below:
    B.profiles[b,] <- mwe.k.profile(data[all.idx[-seq(fold.idx[b,1],fold.idx[b,2],1)],],
                                    data.test,
                                    k.max)$k.profile[1:(N-fold.len)] / nrow(data.test)
  }

  # note: as long as all folds are of equal size it does not matter whether:
  #   a)  we first compute the "% of observations outside" value per fold and then average over
  #       folds or
  #   b)  first sum up all "counts of obs outside" and the divide by total N
  B.profile <- apply(B.profiles, 2, mean, na.rm=T) # option a), (1 - coverage): % of rows outside conf band
  #B.profile <- 1 - colSums(B.profiles, na.rm=T)/N # option b), (1 - coverage): % of rows outside conf band

  #   plot(B.profile, ylim=c(0,1), type="l")
  #   lines(matrix(c(0,0,245,1),nrow=2,byrow=T))
  #   lines(matrix(c(0,alpha,245,alpha),nrow=2,byrow=T))
  #   browser()
  return(list(profile = 1 - B.profile, #report coverage of the band
              profile.mat = B.profiles,
              profile.k = seq(0,length(B.profile)-1,1),
              fold.idx = fold.idx))
}


#' Test set "k profile" using Minimum Width Envelope greedy algorithm
#'
#' @description
#' Determines the iterations k at which observations in the test set become extreme
#' i.e. do not belong to the confidence band anymore. The result is an (1,N) vector
#' of integers {0,1,...,nrow(data.test} identifying how many of the test set observations
#' are extreme at a given k=[1,...,N].
#'
#' This funcition is needed in the mwe.fold() algorithm.
#'
#' @param data.rs (N,M) numeric matrix, N samples, M dimensions/variables. The samples define the
#'                   empirical cdfs of the m variables.
#' @param data.test (N.test,M) numeric matrix, Test dataset
#' @param k.max (1,1) integer, Threshold for ending the computation. When more than 'k.max'
#'                   of the observations in 'data.test' are outside confidence band the profile
#'                   computation ends. Can be used to speed up computations when only the beginning of the profile
#'                   curve is of interest.
#'
#' @return A list with elements:
#' \item{k.profile}{(1,N) vector of integers [0,1,...,nrow(data.test)]}
#'
#' @export
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

        min.df <- min_n_set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
        #Note: rmng.row.arr contains the previous 2nd smallest row -> hence look for 2nd smallest el
        CER$n.min[2,m] = rmng.row.arr[min.df$ind[2]]
        CER$r.min[2,m] = min.df$value[2]

      } else if (CER$n.min[2,m] == n.opt){
        min.df <- min_n_set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
        CER$n.min[2,m] = rmng.row.arr[min.df$ind[2]]
        CER$r.min[2,m] = min.df$value[2]

      } else if (CER$n.max[1,m] == n.opt){
        CER$n.max[1,m] = CER$n.max[2,m]
        CER$r.max[1,m] = CER$r.max[2,m]
        CER$env.max[m] = data.rs[CER$n.max[1,m],m]

        max.df <- max_n_set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
        CER$n.max[2,m] = rmng.row.arr[max.df$ind[2]]
        CER$r.max[2,m] = max.df$value[2]

      } else if (CER$n.max[2,m] == n.opt){
        max.df <- max_n_set(data.rank[rmng.row.arr,m],k=2) #O(N.mwe)
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
