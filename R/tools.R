
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


## repmat (as in Matlab)
# http://waxworksmath.com/Authors/G_M/Hastie/Code/Chapter4/repmat.R
repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}
