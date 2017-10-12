#' Other multivariate confidence intervals and related tools

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
#'
#' @export
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
#'
#' @export
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
#'
#' @export
cb_quantile <- function(dmat, k, quantile.fun = min_width_quantile){
  start_time <- Sys.time()

  cb <- apply(dmat, 2, quantile.fun, k = k)

  list( cb = cb,
        k = k,
        quantile.fun = quantile.fun,
        wall_time_taken = as.numeric(difftime(Sys.time(), start_time, units = "sec")) )
}


##########################################################################################
## Misc
##########################################################################################

#' Data envelope and related statistics
#'
#' @description
#'
#' @param data  [N,M] numeric matrix, N samples, M variables/dimensions.
#'
#' @return A list with elements:
#' \item{envelope}{[2,M] matrix, lower (row 1) and upper (row 2) envelopes}
#' \item{sum.width}{[1,1] numeric, Envelope width (sum of widths over the M dimensions)}
#' \item{mean.width}{[1,1] numeric, Mean envelope width (mean of widths over the M dimensions)}
#' \item{min/max}{you get it...}
#'
#' @export
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


#' Find observations that are outside confidence intervals in at least one variable
#'
#' @description
#'
#' @param data.train (n,m) numeric matrix, n samples, m variables. The samples define the empirical cdfs of the m variables. Variables correspond to multiple hypotheses.
#' @param cb.low (1,m) numeric vector, lower CI limit for each of the dimensions
#' @param cb.high (1,m) numeric vector, upper CI limit for each of the dimensions
#'
#' @return [M,1] logical, True if the observation (row) exits the mvci
#'
#' @export
row.outside.ci <- function(data.train, cb.low, cb.high){
  N = nrow(data.train)
  data.out.match <- (data.train < repmat(matrix(cb.low, nrow=1, byrow=T), N, 1)) |
    (repmat(matrix(cb.high, nrow=1, byrow=T), N, 1) < data.train)
  return(0 < rowSums(data.out.match))
}


#' #' Count of cb violations in vector x
#' cb.violations <- function( x, cb.low, cb.high ) {
#'   sum( sapply( 1:length(cb.high), function(i) cb.low[i] > x[i] || cb.high[i] < x[i] ) )
#' }
#'
#' #' Does vector x violate cb with criterion L
#' cb.contains <- function( x, cb.low, cb.high, L ) {
#'   cb.violations( x, cb.low, cb.high ) <= L
#' }
#'
#' #' Total number of rows in X that do not violate cb
#' total.contained <- function( X, cb.low, cb.high, L ) {
#'   sum( apply( X, 1, function(x) cb.contains( x, cb.low, cb.high, L ) ) )
#' }


#' Compute Pr[V(x|x_u,x_l)<=L] using X_train, K, L and X_test
#'
#' @description
#' TODO: L or l, cb or mvci, K or k, etc.
#'
#' @param cbfun function, Function to compute mvci with, should take only X_train as input
#' @param X_train [N0,M] numeric matrix, Training data, N samples, M variables.
#' @param X_test [N1,M] numeric matrix, Test data, N samples, M variables.
#' @param L [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#'
#' @return [1,1] numeric, percentage of observations in X_test that are within the mvci produced by X_train
#'
#' @export
prc.inside.cb <- function(cbfun, X_train, X_test, L){

  # debug:
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


#' Width of the mvci
#'
#' @description
#'
#' @param cb list, mvci specifications produced by one of the mvci functions
#'
#' @return [1,1] numeric, Width of the mvci (the MWE target cost criterion)
#'
#' @export
cb_area <- function( cb ) {
  ciup   <- cb$data[cb$upmask]
  cidown <- cb$data[cb$downmask]

  sum( ciup - cidown )
}
