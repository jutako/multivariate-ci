
#' k largest values of x (an extension of base::max())
#'
#' @description
#' Searches for maximum values up to k:th largest value. See k.min() for details.
#'
#' @param x (1,N) numeric, Data vector
#' @param k (1,1) integer, Up to how far to search
#'
#' @return [1,J] A data.frame with information on the k largest values of x. See k.min() for details.
#'
#' @export
max_n_set <- function(x,k){
  if (k <= 0){ stop("k.min: too small k") }
  if (length(x) < k){ stop("k.min: too large k") }
  res <- min_n_set(-x,k)
  res$value <- -res$value
  return(res)
}


#' Find n:th largest value of x
#'
#' @description
#' Find n:th largest value of x
#' From: http://stackoverflow.com/questions/2453326/fastest-way-to-find-second-third-highest-lowest-value-in-vector-or-column
#'
#' @param x [1,N] numeric, Data vector
#' @param n [1,1] integer, Up to how far to search
#'
#' @return [1,1] numeric, n:th largest value of x
#'
#' @export
max_n <- function(x, n){
  len <- length(x)
  if(n > len){
    warning('N greater than length(x).  Setting N=length(x)')
    n <- length(x)
  }
  sort(x, partial = len-n+1)[len-n+1]
}


#' k smallest values of x (an extension of base::min())
#'
#' @description
#' Searches for minimum values up to k:th smallest value in O(k*N) time.
#'
#' @param x (1,N) numeric, Data vector
#' @param k (1,1) integer, Up to how far to search
#'
#' @return A data.frame with k rows and columns:
#' \item{k}{1 for min, 2 for 2nd smallest, etc. ...}
#' \item{ind}{position within x}
#' \item{value}{x[ind]}
#' so that x[res$ind] = res$value are the k smallest values, starting with the smallest.
#'
#' @export
min_n_set <- function(x,k){
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


#' Find n:th smallest value of x
#'
#' @description
#' Find n:th smallest value of x
#'
#' @param x [1,N] numeric, Data vector
#' @param n [1,1] integer, Up to how far to search
#'
#' @return [1,1] numeric, n:th smallest value of x
#'
#' @export
min_n <- function(x, n){
  -max_n(-x, n)
}


#' Create a large matrix from copies of a small one
#'
#' @description
#' Functions like repmat() in Matlab:
#' http://waxworksmath.com/Authors/G_M/Hastie/Code/Chapter4/repmat.R
#'
#' @param X [a,b] numeric, A matrix
#' @param i [1,1] integer, Number of repeats over rows
#' @param k [1,1] integer, Number of repeats over cols
#'
#' @return [a*n, b*m] numeric, A matrix created by repeating X
#'
#' @examples
#' repmat(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE), 4, 3)
#'
#' @export
repmat = function(X,i,k){
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*k)), mx*i, nx*k, byrow = TRUE)
}
