# Code that implements the minimum interval quantile MVCI (multivariate confidence interval)

#' Computes minimum width quantile CI separately for each column of matrix X i.e.
#' separately for each variable
#' 
#' Adjusts the number of observations to exclude (i.e. the percentage a in function
#' f()) such that the total number of rows insice the confidence area mathces n-k
mvci_min_width_quantile <- function( X, k, l ) {
  X <- as.matrix( X )
  n <- nrow(X)
  f <- function(a) {
    ci <- get_ci( X, a )
    total_contained( ci, X, l ) - (n-k)
  }
  aopt <- uniroot( f, c(0,0.5) )$root
  get_ci( X, aopt )
}


# A wrapper for backward combatibility
mvci_uniform_percentiles <- mvci_min_width_quantile

#' Mininum width quantile (1-a)*100% CI for a data matrix
#' 
#' Applies smallest_ci() to each column of X
get_ci <- function( X, a ) {
  k <- floor( nrow(X)*a )
  foo <- apply( X, 2, function(x) smallest_ci( x, k ) )
  list( xu=foo[2,], xl=foo[1,] )
}


#' Mininum width quantile CI for a vector of data
#' Returns (1 - k/length(x)) * 100% CI for a data vector using minimum width quantiles
#' 
#' todo: same as min_width_quantile()
smallest_ci <- function(x,k) {
  n <- length(x)
  x <- sort(x)
  i <- which.min(x[(n-k):n]-x[1:(k+1)])
  c(x[i],x[n+i-k-1])
}


#' Returns the number of points in which x lies outside ci
ci_violations <- function( ci, x ) {
  sum( ci$xu - x < 0 ) + sum( x - ci$xl < 0 )
}


#' Returns TRUE if ci contains x with respect to relaxation l and FALSE otherwise
ci_contains <- function( ci, x, l ) {
  ci_violations( ci, x ) <= l
}


#' Returns the number of rows in matrix X that are within ci with respect
#' to relaxation l
total_contained <- function( ci, X, l ) {
  sum( apply( X, 1, function(x) ci_contains( ci, x, l ) ) )
}
